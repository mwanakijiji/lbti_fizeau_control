# This tries out the distribution of processes across multiple cores, and
# writing the results without collisions to a single file

from modules import *
from modules import plot_injected_retrieved
import multiprocessing as mp
import numpy as np
import pandas as pd
import pyfits
import glob
import os
import time

fn = './results_trial1_opd_tip_tilt.txt'

def worker(file_name, q):
    # reads in a single FITS file, calculates some quantities, and returns them

    string_basename = os.path.basename(file_name)

    f = pyfits.open(file_name)
    header = f[0].header # contains injected values

    img_before_padding_before_FT = np.copy(f[0].data)
    print(type(img_before_padding_before_FT))

    size_y = np.shape(img_before_padding_before_FT)[0]
    size_x = np.shape(img_before_padding_before_FT)[1]

    # locate PSF (note this is only to precision of a pixel)
    psf_loc = find_grism_psf(img_before_padding_before_FT, sig=5, length_y=5)
    x_shift_retrieve = np.subtract(psf_loc[1],0.5*np.shape(img_before_padding_before_FT)[1])
    y_shift_retrieve = np.subtract(psf_loc[0],0.5*np.shape(img_before_padding_before_FT)[0])

    # take FFT
    amp, arg = fft_img(img_before_padding_before_FT).fft(padding=200, mask_thresh=0)

    # analyze FFTs
    fftInfo_amp = fftMask(amp,sci_wavel,plateScale,fyi_string=" FFT amp")
    fftInfo_arg = fftMask(arg,sci_wavel,plateScale,fyi_string=" FFT phase")

    # thresholds
    fft_ampl_high_freq_lowlimit = 2.4e5 # min threshold for good fringe visibility
    fft_ampl_low_freq_lowlimit = 1.4e6 # min threshold for acceptable AO correction
    fft_phase_vec_high_freq_highlimit = 5 # max threshold for Airy overlap
    std_lowFreqPerfect_lowlimit = 10 # max threshold for Airy overlap
    phase_normVec_highFreqPerfect_R_x = 100 # max threshold for phase of high-freq fringes

    # find needed correction to FPC PL setpoint (degrees in K-band)
    #sci_to_K = np.divide(sci_wavel,2.2e-6) # factor to convert degrees in sci to degrees in K
    #corrxn_pl = -fftInfo_arg["med_highFreqPerfect_R"].values[0]*(180./np.pi)*sci_to_K
    #alpha_high_freq = [x_grad_perf_high_R, y_grad_perf_high_R] # gradient high freq lobe of PTF in x and y: [a,b]

    # state size of FFTed image
    Ny = np.shape(amp.data)[0]
    Nx = np.shape(amp.data)[1]
    # state gradient of PTF slope in x and y
    x_grad_perf_high_R = np.median(fftInfo_arg["normVec_highFreqPerfect_R_x"])
    y_grad_perf_high_R = np.median(fftInfo_arg["normVec_highFreqPerfect_R_y"])
    x_grad_perf_lowfreq = np.median(fftInfo_arg["normVec_lowFreqPerfect_x"])
    y_grad_perf_lowfreq = np.median(fftInfo_arg["normVec_lowFreqPerfect_y"])
    alpha_high_freq = [x_grad_perf_high_R, y_grad_perf_high_R] # gradient high freq lobe of PTF in x and y: [a,b]
    alpha_low_freq = [x_grad_perf_lowfreq, y_grad_perf_lowfreq] # same, in low freq lobe
    alpha_mean = np.mean([alpha_high_freq,alpha_low_freq],axis=0) # tip-tilt corrections should be based on gradients common to lobes (see Spalding+ SPIE 2018, Table 3)
    corrxn_tt = needed_tt_setpt_corrxn(alpha=alpha_mean,PS=plateScale_LMIR,Nx=Nx,Ny=Ny) # (x,y)

    opd_retrieve = fftInfo_arg["med_highFreqPerfect_R"]*((180./np.pi)/360.)*sci_wavel*1e6
    tip_retrieve = -corrxn_tt[1] # y (negative, because I want to know the retrieved value, not the correction)
    tilt_retrieve = -corrxn_tt[0] # x (ditto)

    opd_inject = header["OPD_UM"]
    tip_inject = header["TIPY_MAS"]
    tilt_inject = header["TILTXMAS"]

    # x,y position of PSF may or may not be in the header
    if ("X_SHIFT" in header.keys()):
        x_shift_inject = header["X_SHIFT"]
        y_shift_inject = header["Y_SHIFT"]
    else:
        x_shift_inject = -999
        y_shift_inject = -999

    string_opd_retrieve = str(opd_retrieve)
    string_tip_retrieve = str(tip_retrieve)
    string_tilt_retrieve = str(tilt_retrieve)
    string_x_shift_retrieve = str(x_shift_retrieve)
    string_y_shift_retrieve = str(y_shift_retrieve)

    string_opd_inject = str(opd_inject)
    string_tip_inject = str(tip_inject)
    string_tilt_inject = str(tilt_inject)
    string_x_shift_inject = str(x_shift_inject)
    string_y_shift_inject = str(y_shift_inject)

    print("opd retrieved: " + str(opd_retrieve) + "; injected: " +str(opd_inject))
    print("tip retrieved: " + str(tip_retrieve) + "; injected: " +str(tip_inject))
    print("tilt retrieved: " + str(tilt_retrieve) + "; injected: " +str(tilt_inject))

    '''
    # High-freq fringe visibility (median)
    print("--------------------------")
    print("Median of ampl of high freq lobe:")
    print(fftInfo_amp["med_highFreqPerfect_R"])
    # TO CORRECT: MOVE THE SPC_TRANS TO FIND CENTER OF COHERENCE ENVELOPE

    # High-freq phase gradient
    print("--------------------------")
    print("Median of phase (science PSF):")
    print(fftInfo_arg["med_highFreqPerfect_R"]*(180./np.pi))
    print("Pathlength correction needed (um):")
    print(fftInfo_arg["med_highFreqPerfect_R"]*((180./np.pi)/360.)*sci_wavel*1e6)
    print("Phase gradient in x of high freq in PTF:")
    print(fftInfo_arg["normVec_highFreqPerfect_R_x"])
    print("Phase gradient in y of high freq in PTF:")
    print(fftInfo_arg["normVec_highFreqPerfect_R_y"])
    print("--------------------------")
    '''

    # save images to check
    '''
    hdu = pyfits.PrimaryHDU(img_before_padding_before_FT)
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("junk_seen_prepradding_" + string_basename, clobber=True)
    hdu = pyfits.PrimaryHDU(amp.data)
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("junk_ampPE_" + string_basename + ".fits", clobber=True)
    hdu = pyfits.PrimaryHDU(arg.data)
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("junk_argPE_" + string_basename + ".fits", clobber=True)
    # save fyi FITS files to see the masks, etc.
    hdu = pyfits.PrimaryHDU(fftInfo_amp["sciImg1"])
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("log_images/fft_amp_masked_region_1_" + string_basename + ".fits", clobber=True)
    hdu = pyfits.PrimaryHDU(fftInfo_amp["sciImg2"])
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("log_images/fft_amp_masked_region_2_" + string_basename + ".fits", clobber=True)
    hdu = pyfits.PrimaryHDU(fftInfo_amp["sciImg3"])
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("log_images/fft_amp_masked_region_3_" + string_basename + ".fits", clobber=True)
    hdu = pyfits.PrimaryHDU(fftInfo_amp["sciImg4"])
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("log_images/fft_amp_masked_region_4_" + string_basename + ".fits", clobber=True)
    hdu = pyfits.PrimaryHDU(fftInfo_arg["sciImg1"])
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("log_images/fft_arg_masked_region_1_" + string_basename + ".fits", clobber=True)
    hdu = pyfits.PrimaryHDU(fftInfo_arg["sciImg2"])
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("log_images/fft_arg_masked_region_2_" + string_basename + ".fits", clobber=True)
    hdu = pyfits.PrimaryHDU(fftInfo_arg["sciImg3"])
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("log_images/fft_arg_masked_region_3_" + string_basename + ".fits", clobber=True)
    hdu = pyfits.PrimaryHDU(fftInfo_arg["sciImg4"])
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("log_images/fft_arg_masked_region_4_" + string_basename + ".fits", clobber=True)
    '''

    with open(fn, 'rb') as f:
        size = len(f.read())
    list1 = string_basename, string_opd_inject, string_tip_inject, string_tilt_inject, string_x_shift_inject, string_y_shift_inject, \
            string_opd_retrieve, string_tip_retrieve, string_tilt_retrieve, string_x_shift_retrieve, string_y_shift_retrieve
    res = ','.join(list1) # strip extra parentheses etc.
    q.put(res)

    return res

def listener(q):
    '''listens for messages on the q, writes to file. '''

    with open(fn, 'w') as f:
        while 1:
            m = q.get()
            if m == 'kill':
                break
            f.write(str(m) + '\n')
            f.flush()

def main():
    '''
    Grab frames from a directory and retrieve the aberrations in the PSF
    '''

    # choose the directory
    path_stem = "./synthetic_fizeau/trial1_opd_tip_tilt/"
    #path_stem = "./synthetic_fizeau/trial2_opd/"
    #path_stem = "./synthetic_fizeau/trial3_tip/"
    #path_stem = "./synthetic_fizeau/trial4_tilt/"
    #path_stem = "./synthetic_fizeau/trial5_yx/"
    #path_stem = "./synthetic_fizeau/trial6_opd_tip_tilt_yx/"
    #path_stem = "/vol_c/synthetic_fizeau/trial2_opd/"

    # get list of file names together
    files_list = glob.glob(path_stem + "*.fits")
    print("Files to analyze: ")
    print(files_list)

    #must use Manager queue here, or will not work
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() + 2)

    #put listener to work first
    watcher = pool.apply_async(listener, (q,))

    #fire off workers
    jobs = []
    for i in range(len(files_list)):
        file_name_fits = files_list[i]
        job = pool.apply_async(worker, (file_name_fits, q))
        jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    #now we are done, kill the listener
    q.put('kill')
    pool.close()
    pool.join()


if __name__ == "__main__":
    main() # analysis of frames
    #plot_injected_retrieved.plot_analysis(csv_file = fn) # plotting of analysis
