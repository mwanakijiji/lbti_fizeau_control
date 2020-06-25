# This tries out the distribution of processes across multiple cores, and
# writing the results without collisions to a single file

from modules import *
import multiprocessing as mp
import numpy as np
import pyfits
import glob
import os
import time

fn = './junk.txt'

def worker(file_name, q):
    # reads in a single FITS file, calculates some quantities, and returns them

    string_basename = os.path.basename(file_name)

    f = pyfits.open(file_name)
    img_before_padding_before_FT = np.copy(f[0].data)
    print(type(img_before_padding_before_FT))

    size_y = np.shape(img_before_padding_before_FT)[0]
    size_x = np.shape(img_before_padding_before_FT)[1]

    # locate PSF
    psf_loc = find_grism_psf(img_before_padding_before_FT, sig=5, length_y=5)
    print(psf_loc)

    # take FFT
    amp, arg = fft_img(img_before_padding_before_FT).fft(padding=200, mask_thresh=0)

    # analyze FFTs
    fftInfo_amp = fftMask(amp,sci_wavel,plateScale,fyi_string=" FFT amp")
    fftInfo_arg = fftMask(arg,sci_wavel,plateScale,fyi_string=" FFT phase")

    hdu = pyfits.PrimaryHDU(fftInfo_amp["sciImg1"])
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("fft_amp_masked_region_1_" + string_basename + ".fits", clobber=True)


    string_size_y = str(size_y)
    string_size_x = str(size_x)

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
    list1 = string_basename, string_size_y, string_size_x
    res = ', '.join(list1) # strip extra parentheses etc.
    q.put(res)

    return res

def listener(q):
    '''listens for messages on the q, writes to file. '''

    with open(fn, 'w') as f:
        while 1:
            m = q.get()
            if m == 'kill':
                f.write('killed')
                break
            f.write(str(m) + '\n')
            f.flush()

def main():
    '''
    Grab frames from a directory and retrieve the aberrations in the PSF
    '''

    path_stem = "./synthetic_fizeau/trial1_opd_tip_tilt/"

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
   main()
