import sys, os, string, time, pdb, copy, pyfits
import numpy as np
#import np.ma as ma
import pandas as pd
from pyindi import *
import scipy
from scipy import ndimage, sqrt, stats, misc, signal
import pyfits
import pickle
import math
import glob
from lmircam_tools import *
from lmircam_tools import process_readout

def live_opd_correction_fizeau_grism(integ_time, mode = "science"):
    ''' 
    Measures the angle of fringes in grism mode, and calculates
    and implements the needed SPC translation stage movement

    INPUTS:
    mode: testing or science
    '''

    # set some approximate parameters of the observed grism PSF
    sig = 5 # sigma of Gaussian profile in x (in pix)
    length_y = 200 # length in y of the psf (in pix)

    counter_num = 0 # for counting number of analyzed PSFs
    angle_val_array = [] # initialize array which will store values from individual frames

    # read in any new images written out to a directory
    files_start = glob.glob(dir_to_monitor + "*.fits") # starting list of files
    num_psfs_to_analyze = 10 # number of PSFs to sample

    while counter_num < num_psfs_to_analyze:

        # if we're just reading in fake FITS files, monitor a directory
        # and read in a new file
        if (mode == "fake_fits"):

            while True:

                time_start = time.time()
                time.sleep(del_t)

                # filename for pickled FFT info
                fft_pickle_write_name = "pickled_info/fft_info_"+str("{:0>2d}".format(counter_num))+".pkl"

                # check to see if there were new files from last check
                files_later = glob.glob(dir_to_monitor + "/*.fits")

                # are there new files?
                new_list = np.setdiff1d(files_later,files_start)

                time_start = time.time() # start timer

                # if there are no new files, cycle back through
                if (len(new_list) <= 2):
                    continue

                # if there are new files, read one of them in
                elif (len(new_list) > 2):

                    # reassign these files to be next starting point
                    files_start = files_later

                    # filename of second-newest file (in case the newest is still being written)
                    second_newest = sorted(files_later)[-2]

                    # read in file
                    f = pyfits.open(second_newest)
                    file_name_full = os.path.basename(second_newest) # filename without the path
                    file_name_base = os.path.splitext(file_name_full)[0] # filename without the extension, too

                    # break out of the 'while True'
                    break

        # if we're reading in real frames from LMIR, catch them as they are read out
        elif ((mode == "az_source") or (mode == "science")):

            f = pi_fiz.getFITS("fizeau.roi_image.file", timeout=60)

            # FFT mag image (this has to wait for a new frame)
            fftw_amp = pi_fiz.getFITS("fizeau.mag_image.file", timeout=60)

            # FFT phase image (this has to wait for a new frame(
            fftw_phase = pi_fiz.getFITS("fizeau.phase_image.file", timeout=60)

        # get the right image slice
        if (mode == "fake_fits"):
            if (np.ndim(f[0].data) > 2):
                image = f[0].data[-1,:,:] # images from LMIRcam (> summer 2018) are cubes of nondestructive reads
            else:
                image = np.squeeze(f[0].data)
        elif ((mode == "az_source") or (mode == "science")):
            image = np.copy(f[0].data)

        # if this is a fake fits file, do a quick-and-dirty background subtraction
        if (mode == "fake_fits"):
            image = process_readout.processImg(image,"median")

        # save detector image to check (overwrites previous)
        #hdu = pyfits.PrimaryHDU(image)
        #hdulist = pyfits.HDUList([hdu])
        #hdu.writeto("junk_other_tests/junk_test_image_seen.fits", clobber=True)

        # determine grism Fizeau PSF center
        center_grism = find_grism_psf(image, sig, length_y) # locate the grism PSF center (THIS IS IN OVERLAP_PSFS.PY; SHOULD IT BE IN INIT?)
        # if we are reading in fake FITS files, we may have to just set the location
        if (mode == "fake_fits"):
            center_grism = psf_loc_fake
            # cut out the grism image (best to have rectangle, rather than square cutout)
            #img_before_padding_before_FT = np.copy(image)
            img_before_padding_before_FT = image[center_grism[0]-int(200):center_grism[0]+int(200),center_grism[1]-int(100):center_grism[1]+int(100)]
        elif ((mode == "az_source") or (mode == "science")):
            img_before_padding_before_FT = np.copy(f[0].data)
            file_name_base = str(int(time.time())) + "_" + str(counter_num) # basename has already been defined for fake_fits

        # take FFT; no padding for now
        ## ## DO I WANT A SMALL CUTOUT OR THE ORIGINAL IMAGE?
        AmpPE, ArgPE = fft_img(img_before_padding_before_FT).fft(padding=0)

        # this is a kludge for slipping in the INDI FFT amplitude in place of the Python one
        # (the phase has a checkerboard pattern until Paul fixes it, so Im just going to keep 
        # the Python amplitude)
        if ((mode == "az_source") or (mode == "science")):
            AmpPE = fftw_amp[0]
        # save image to check
        hdu = pyfits.PrimaryHDU(img_before_padding_before_FT)
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/img_seen_prepradding_" + file_name_base + ".fits", clobber=True)
        print('4')

        # test: see what the FFT looks like
        #if (mode == "fake_fits"):
        #    hdu = pyfits.PrimaryHDU(AmpPE.data)
        #else: # if we are using FFTW FFT magnitude
        #    hdu = pyfits.PrimaryHDU(fftw_amp[0].data)
        hdu = pyfits.PrimaryHDU(AmpPE.data)
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/fft_amp_" + file_name_base + ".fits", clobber=True)
        print('4b')
        hdu = pyfits.PrimaryHDU(ArgPE.data)
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/fft_arg_" + file_name_base + ".fits", clobber=True)
        print('5')
        # find angle of fringes
        # is there an off-center dot in FFT amplitude?
        # blot out low-frequency center of FFT ampl
        center_masked_data = AmpPE.data
        center_masked_data[int(0.5*np.shape(center_masked_data)[0])-20:int(0.5*np.shape(center_masked_data)[0])+20,
                           int(0.5*np.shape(center_masked_data)[1])-20:int(0.5*np.shape(center_masked_data)[1])+20] = np.nan
        dot_loc = find_airy_psf(center_masked_data)
        #dot_loc = find_grism_psf(np.multiply(amp.data,amp.mask), 5, 5)

        print('6')
        # save image to check
        hdu = pyfits.PrimaryHDU(center_masked_data)
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/img_when_centroiding_fringe_angle_" + file_name_base + ".fits", clobber=True)

        print('7')
        print("Analyzing "+file_name_base)
        print("Dot location (y,x):")
        print(dot_loc)
        print("Dot angle (deg CCW from +x):")
        y_comp = dot_loc[0]-0.5*np.shape(center_masked_data)[0]
        x_comp = dot_loc[1]-0.5*np.shape(center_masked_data)[1]
        angle_val = math.atan2(y_comp,x_comp)*180./np.pi
        print(angle_val)

        # concatenate the angle to the array of single-frame samples
        angle_val_array = np.concatenate((angle_val_array, [angle_val]))

        # increase the counter
        counter_num += 1
        print("-----------------")

    #############################################################################
    # IMPLEMENTATION

    # find median angle
    # ... and first remove zeros from the array (which indicates the fringes weren't sensed)
    angle_val = np.nanmedian(angle_val_array[np.nonzero(angle_val_array)])
    print("Median angle value is " + str(angle_val))

    # as found by using OPD scans in grism mode in 2018A and 2018B, it appears that, for the Lgrism6AR,
    # for every +degree in the CW direction that the FFT amplitude high-freq node is,
    # need to move the SPC_Trans in the NEGATIVE direction by 144 counts
    movement_per_pos_degree = -144
    diff_movement_total_cts = angle_val * movement_per_pos_degree # units of counts
    diff_movement_total_opd = 2. * np.divide(diff_movement_total_cts,50.) # factor of 2: linear to OPD; factor of 1/50: counts to um

    # correct with the SPC translation stage: Ubcs.SPC_Trans.command=>N
    # note factor of 10; command is in relative movement of 0.1 um
    print("----------------------------------------------------------------")
    if ((mode == "az_source") or (mode == "science")):
        print("Moving SPC_Trans for large OPD movement of "+str(int(diff_movement_total_opd))+" um or "+str(diff_movement_total_cts)+" translation counts")
        spc_trans_position_pre = pi.getINDI("Ubcs.SPC_Trans_status.PosNum") # translation stage before command (absolute position, 0.02 um)
        spc_trans_position_now = spc_trans_position_pre
        spc_trans_position_goal = spc_trans_position_pre + diff_movement_total_cts
        pi.setINDI("Ubcs.SPC_Trans.command=>"+'{0:.1f}'.format(diff_movement_total_cts))
        time.sleep(5.0) # this pause may be necessary to avoid confusing INDI; SPC very sticky 2019 05 19, E.S.
        # wait for the SPC_Trans movement to finish
        while (spc_trans_position_now != spc_trans_position_goal):
            spc_trans_position_now = pi.getINDI("Ubcs.SPC_Trans_status.PosNum") # check status of SPC_Trans
            time.sleep(0.5) # wait a moment
    else:
        print("Not moving SPC_Trans, since this is in testing mode. But calculated needed movement is " + str(int(diff_movement_total_opd))+" um or "+str(diff_movement_total_cts)+" translation counts")

    return



def find_optimal_opd_fizeau_grism(integ_time, mode = "science"):
    #Takes well-overlapped grism PSFs and dials the optical path
    #difference such that barber-pole fringes become vertical

    #INPUTS:
    #psf_location: location of the Fizeau grism PSF
    #mode: testing or science

    # set some approximate parameters of the observed grism PSF
    sig = 5 # sigma of Gaussian profile in x (in pix)
    length_y = 200 # length in y of the psf (in pix)

    raw_input("User: remove the Blank in FW4, then press return when done")

    # initialize dataframe for pathlength and residual data
    df = pd.DataFrame(columns=['step',\
				'spc_trans_position',\
				'spc_trans_position_opd_um',\
				'hpc_piezo_piston',\
				'hpc_piezo_piston_opd_um',\
				'fpc_piezo_piston',\
				'fpc_piezo_piston_opd_um',\
				'resid'])
    # loop over OPD steps (unitless here) around where we start
    step_start = -1
    step_stop = 4
    for opd_step in range(step_start,step_stop):

	step_size_opd = 10. # step size per opd_step count (um, total OPD)

        if (mode != "fake_fits"):
            print("Taking a background-subtracted frame")

            # ROI image
            f = pi_fiz.getFITS("fizeau.roi_image.file", timeout=60)

            # FFT mag image (this has to wait for a new frame)
            fftw_amp = pi_fiz.getFITS("fizeau.mag_image.file", timeout=60)

            # FFT phase image (this has to wait for a new frame(
            fftw_phase = pi_fiz.getFITS("fizeau.phase_image.file", timeout=60)

	elif (mode == "fake_fits"):
            ## ## CHANGE THIS TO MONITOR A DIRECTORY
            f = pyfits.open("test_fits_files/test_frame_grismFiz_small.fits")

        if (np.ndim(f[0].data) > 2):
            imgb4 = f[0].data[-1,:,:] # images from LMIRcam (> summer 2018) are cubes of nondestructive reads
        else:
            imgb4 = np.squeeze(f[0].data)

        if (mode == "fake_fits"):
            # return background-subtracted, bad-pix-corrected image
            image = process_readout.processImg(imgb4, 'median', background = True)
        elif ((mode == "az_source") or (mode == "science")):
            # just return bad-pix-corrected image
            image = process_readout.processImg(imgb4, 'median', background = False)


        # determine grism Fizeau PSF center
        center_grism = find_grism_psf(image, sig, length_y) # locate the grism PSF center (THIS IS IN OVERLAP_PSFS.PY; SHOULD IT BE IN INIT?)

        # cut out the grism image
        img_before_padding_before_FT = np.copy(image)
        ''' 
        img_before_padding_before_FT = image[center_grism[0]-int(0.5*length_y):center_grism[0]+int(0.5*length_y),
                                     center_grism[1]-2*sig:center_grism[1]+2*sig]
        '''

        # take FFT; no padding for now
        ## ## DO I WANT A SMALL CUTOUT OR THE ORIGINAL IMAGE?
	AmpPE, ArgPE = fft_img(img_before_padding_before_FT).fft(padding=0)

        # this is a kludge for slipping in the INDI FFT amplitude (the phase has a checkerboard pattern until Paul fixes it) in place of the Python one
        if ((mode == "az_source") or (mode == "science")):
            AmpPE = ma.masked_where(fftw_amp == np.nan, fftw_amp, copy=False)

	# save fyi FITS files
        '''
	hdu = pyfits.PrimaryHDU(image)
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto('junk_test_original_grism.fits', clobber=True)
    	hdu = pyfits.PrimaryHDU(AmpPE.data)
    	hdulist = pyfits.HDUList([hdu])
    	hdu.writeto('junk_test_amp_grism.fits', clobber=True)
        if (mode == "fake_fits"):
    	    hdu = pyfits.PrimaryHDU(ArgPE.data)
        else: # if we are using FFTW FFT magnitude
            hdu = pyfits.PrimaryHDU(fftw_phase)
    	hdulist = pyfits.HDUList([hdu])
    	hdu.writeto('junk_test_arg_grism.fits', clobber=True)
	'''

        #############################################################
        ## ## BEGIN OPTIMIZATION OF PATHLENGTH BASED ON GRISM FFTS

	# LATER DO THE FOLLOWING; FOR NOW, JUST DIAL IN BOTH DIRECTIONS AND FIT A POLYNOMIAL
        # if direction of fringes is apparent --> dial OPD until past the point
        # where they are vertical, then fit a polynomial
        # if no direction of fringes is apparent --> dial one way 10 microns, then
	# jump back and go the other way 10 microns

	# The idea here is to
	# 1. Integrate over one side (left or right) of the 2D FFT of the grism image
	# 2. Take the cross-correlation of the upper and lower part of the above
	# 3. Subtract one side of the correlation from the other, to find asymmetries
	# 4. Insert the median value of the residual to an array
	# 5. Move to the next OPD position and re-take measurements
	# 6. After going through a given range of movement, fit a polynomial to the 
	#     residuals as a function of OPD

	# take right-hand side of FFT, integrate in x
        to_corr = np.sum(AmpPE[:,int(0.5*img_before_padding_before_FT.shape[1]):], axis=1) 
        to_corr_masked = np.copy(to_corr)
	# mask the strong low-frequency power
	to_corr_masked[int(0.5*len(to_corr_masked))-2:int(0.5*len(to_corr_masked))+2] = 0 
	# take cross-correlation
        test_symm = signal.correlate(to_corr_masked, to_corr_masked[::-1], mode='same') 
	# separate into halves
        leftHalf = test_symm[:int(0.5*len(test_symm))]
        rightHalf = test_symm[int(0.5*len(test_symm)):]
        # find residuals between the halves of the cross-correlation
	resid = leftHalf-rightHalf[::-1]

	# add median of residuals, and pathlength (HPC) position, to arrays
	print("----------------------------------------------------------------")
	print("Current OPD step:")
	print(str(opd_step) + ", to end with " + str(np.add(step_start,np.subtract(step_stop,step_start))))

	# conversions
	# SPC_Trans: 1 um trans (2 um OPD) -> 50 counts
	# HPC piezo: 1 um trans (2 um OPD) -> 10 counts
	# FPC piezo: 1 um trans (2 um OPD) -> 10 counts
	spc_trans_position = pi.getINDI("Ubcs.SPC_Trans_status.PosNum") # translation stage (absolute position, 0.02 um)
	hpc_piezo_piston = pi.getINDI("Acromag.HPC_status.Piston") # piezo piston (absolute position, um)
	fpc_piezo_piston = pi.getINDI("Acromag.FPC_status.Piston") # piezo piston (absolute position, um)
	spc_trans_position_opd_um = 2.*np.divide(spc_trans_position,50.) # (OPD, um)
	hpc_piezo_piston_opd_um = 2.*np.divide(hpc_piezo_piston,10.) # (OPD, um)
	fpc_piezo_piston_opd_um = 2.*np.divide(fpc_piezo_piston,10.) # (OPD, um)

	#df.append(pd.DataFrame())
	df_append = pd.DataFrame(np.transpose([[opd_step],\
					[spc_trans_position],\
					[spc_trans_position_opd_um],\
					[hpc_piezo_piston],\
					[hpc_piezo_piston_opd_um],\
					[fpc_piezo_piston],\
					[fpc_piezo_piston_opd_um],\
					[np.median(resid)]]), \
				columns=['step',\
					'spc_trans_position',\
					'spc_trans_position_opd_um',\
					'hpc_piezo_piston',\
					'hpc_piezo_piston_opd_um',\
					'fpc_piezo_piston',\
					'fpc_piezo_piston_opd_um',\
					'resid'])
	df = df.append(df_append, ignore_index=True)
	print(df)
	# now move the HPC to the next step (small steps with piezos)
        # small steps, piezos: Acromag.HPC.Tip=0;Tilt=0;Piston=[step_size_opd];Mode=1
        #hpc_small_step = 0.5*step_size_opd # half the OPD (relative step)
	#hpc_piezo_next_pos = np.add(spc_trans_position, opd_step*hpc_small_step) # piezo command is in absolute position, units of um
	#print("----------------------------------------------------------------")
	#print("Moving HPC for small OPD movement to position "+hpc_piezo_next_pos)
	#pi.setINDI("Acromag.HPC.Tip=0;Tilt=0;Piston="+'{0:.1f}'.format(hpc_piezo_next_pos)+";Mode=1")

	# big steps with the SPC translation stage: Ubcs.SPC_Trans.command=>5
	# note factor of 10; command is in relative movement of 0.1 um
        print("----------------------------------------------------------------")
        if ((mode == "az_source") or (mode == "science")):
	    print("Moving SPC_Trans for large OPD movement of "+str(int(step_size_opd))+" um (translation of "+str(0.5*step_size_opd)+\
		" um, or "+str(50*0.5*step_size_opd)+" counts)")
	    pi.setINDI("Ubcs.SPC_Trans.command=>"+'{0:.1f}'.format(50*0.5*step_size_opd))
        elif (mode == "fake_fits"):
            print("Since this is fake_fits mode, SPC_Trans will not be moved, but it WOULD have been moved for large OPD movement of "+\
		str(int(step_size_opd))+" um (translation of "+str(0.5*step_size_opd)+" um, or "+str(50*0.5*step_size_opd))

    # write data to file for copying and plotting on local machine
    file_name = "resids_test.csv"
    #resids_df = pd.DataFrame(np.transpose([pl_array,resid_array]), columns=["pl","resid"])
    no_return = df.to_csv(file_name)

    # fit a 2nd-order polynomial to the residuals as fcn of SPC position
    if ((mode == "az_source") or (mode == "science")):
    	coeffs = np.polyfit(df["spc_trans_position"], df["resid"], 2)
    elif (mode == "fake_fits"):
        coeffs = [0,0]

    ## WOULD BE COOL TO HAVE A WEBPAGE DISPLAY OF A PLOT OF THE ABOVE
    ## MATPLOTLIB VERSION IS
    ## plt.scatter(pl_array,resid_array)
    ## plt.plot(pl_array, coeffs[2] + coeffs[1]*pl_array+coeffs[0]*pl_array**2)
    ## plt.title('FFT residuals by SPCTRANS position')
    ## plt.xlabel('SPCTRANS')
    ## plt.ylabel('Median of Residuals Between Top and Bottom Halves of FFT')

    # pickle the coefficients and the scan data to be restored by the next function
    # (note this includes the data that was written as a csv for checking)
    with open("pickled_info/scan_objs.pkl", "w") as f:
        pickle.dump([coeffs, df], f)

    return


def implement_optimal_opd(mode = "science"):
    #Take the data from the scan, calculate where OPD=0 is, and implement it by moving the SPC

    #INPUTS:
    #d: dictionary containing
    #    "coeffs": the coefficients of the fit to the residuals
    #    "scan_data": the data on the residuals, and piston positions of moveable beam combiner elements

    # write data to file for copying and plotting on local machine
    ## ## (THIS IS ACTUALLY DEFINED IN FIND_OPTIMAL_FIZEAU_GRISM; THIS FILENAME NEEDS TO BE TRANSFERRED INSTEAD
    file_name = "resids_test.csv"

    # restore the pickle file with the fit coefficients and scan data
    with open("pickled_info/scan_objs.pkl") as f:
        coeffs, df = pickle.load(f)

    # find the minimum; set the HPC path length position accordingly
    y_series = coeffs[2] + coeffs[1]*df["spc_trans_position"]+coeffs[0]*df["spc_trans_position"]**2
    min_y = min(y_series)  # find the minimum y-value
    zero_opd_spc_trans_pos = df.loc[y_series == np.min(y_series)]["spc_trans_position"].values[0] # find the x-value (pathlength) corresponding to min y-value (residuals)

    # is the curve concave up?
    if (coeffs[0] < 0):
        print("Something is wrong-- line fit to residuals is not concave up!!")

    raw_input("Plot "+file_name+" on a local machine to see if OPD of zero has been found. Then proceed with [Enter]. Otherwise, break.")

    # command the SPC_Trans to move to that minimum
    spc_trans_position_now = pi.getINDI("Ubcs.SPC_Trans_status.PosNum") # re-obtain translation stage position
    spc_trans_command = np.subtract(zero_opd_spc_trans_pos,spc_trans_position_now)
    print("spc_trans_position_now:")
    print(spc_trans_position_now)
    if ((mode == "fake_fits") or (mode == "az_source") or (mode == "science")):
        print("spc_trans_command:")
        print(spc_trans_command)
        pi.setINDI("Ubcs.SPC_Trans.command=>"+'{0:.1f}'.format(spc_trans_command))

    # command the HPC piezo to move to that minimum
    #if (mode != "total_passive"):
    #    pi.setINDI("Acromag.HPC.Tip=0;Tilt=0;Piston="+'{0:.1f}'.format(max_x)+";Mode=1")

    # as last step, remove grism
    raw_input("User: remove grism, insert observing filters, and press [Enter]")

    # turn off fizeau flag to avoid problems with other observations
    if (mode != "total_passive"):
        print("De-activating ROI aquisition flag")
        pi_fiz.setINDI("fizeau.enable_run.value=Off")

    return
