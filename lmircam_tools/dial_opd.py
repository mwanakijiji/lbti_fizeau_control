import sys, os, string, time, pdb, copy, pyfits
import numpy as np
#import np.ma as ma
from pyindi import *
import scipy
from scipy import ndimage, sqrt, stats, misc, signal
import pyfits
from lmircam_tools import *
from lmircam_tools import process_readout


def optimize_opd_fizeau_grism(mode = "science"):
    ''' 
    Takes well-overlapped grism PSFs and dials the optical path
    difference such that barber-pole fringes become vertical

    INPUTS: 
    psf_location: location of the Fizeau grism PSF
    mode: "fake_fits": read in fake FITS files (but continue sending LMIR and mirror commands)
          "artif_source": use an artificial source (either laser or pinhole)
          "science": on-sky
    '''

    # set some approximate parameters of the observed grism PSF
    sig = 5 # sigma of Gaussian profile in x (in pix)
    length_y = 200 # length in y of the psf (in pix)

    take_roi_background()

    # loop over OPD steps (unitless here) around where we start
    for opd_step in range(-20,20):

	step_size = 1. # step size per opd_step count (um, total OPD)

	take_roi_background()

        print("Taking a background-subtracted frame")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)

	if (mode == "fake_fits"):
            f = pyfits.open("test_fits_files/test_frame_grismFiz_small.fits")

        imgb4 = f[0].data
        image = process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

        # determine grism Fizeau PSF center
        center_grism = find_grism_psf(image, sig, length_y) # locate the grism PSF center (THIS IS IN OVERLAP_PSFS.PY; SHOULD IT BE IN INIT?)

        # cut out the grism image
        img_before_padding_before_FT = image[center_grism[0]-int(0.5*length_y):center_grism[0]+int(0.5*length_y),
                                     center_grism[1]-2*sig:center_grism[1]+2*sig]

        # take FFT; no padding for now
        ## ## DO I WANT A SMALL CUTOUT OR THE ORIGINAL IMAGE?
	AmpPE, ArgPE = fft_img(img_before_padding_before_FT).fft(padding=0)

	# save fyi FITS files
        hdu = pyfits.PrimaryHDU(image)
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto('junk_test_original_grism.fits', clobber=True)
    	hdu = pyfits.PrimaryHDU(AmpPE.data)
    	hdulist = pyfits.HDUList([hdu])
    	hdu.writeto('junk_test_amp_grism.fits', clobber=True)
    	hdu = pyfits.PrimaryHDU(ArgPE.data)
    	hdulist = pyfits.HDUList([hdu])
    	hdu.writeto('junk_test_arg_grism.fits', clobber=True)

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
        testArray = np.concatenate((testArray,[np.median(resid)]))
	spc_trans_position = pi.getINDI("Ubcs.SPC_Trans_status.PosNum") # translation stage (absolute position, 0.1 um)
        spc_piezo_position = pi.getINDI("Acromag.HPC_status.Piston") # piezo piston (absolute position, um)
	plArray = np.concatenate((plArray,[int(spc_position)])) # pathlength array

	# now move the HPC to the next step (small steps with piezos)
        # small steps, piezos: Acromag.HPC.Tip=0;Tilt=0;Piston=[step_size];Mode=1
        hpc_small_step = 0.5*step_size # half the OPD (relative step)
	hpc_piezo_next_pos = np.add(spc_piezo_position, opd_step*hpc_small_step) # piezo command is in absolute position, units of um
	print("----------------------------------------------------------------")
	print("Moving HPC for small OPD movement to position "+hpc_piezo_next_pos)
	pi.setINDI("Acromag.HPC.Tip=0;Tilt=0;Piston="+'{0:.1f}'.format(hpc_piezo_next_pos)+";Mode=1")

	### IF WE CANT FIND ANYTHING AT ALL, TAKE BIG STEPS WITH THE TRANSLATION STAGE
	# big steps, translation stage: Ubcs.SPC_Trans.command=>5
	# note factor of 10; command is in relative movement of 0.1 um
        #print("----------------------------------------------------------------")
        #print("Moving HPC for large OPD movement of "+str(int(10*0.5*stepSize))+" um")
	## ## pi.setINDI("Ubcs.SPC_Trans.command=>"+'{0:.1f}'.format(10*0.5*stepSize))

    # fit a 2nd-order polynomial to the residuals as fcn of SPC position
    coeffs = np.polyfit(plArray,testArray, 2)

    ## WOULD BE COOL TO HAVE A WEBPAGE DISPLAY OF A PLOT OF THE ABOVE
    ## MATPLOTLIB VERSION IS
    ## plt.scatter(plArray,testArray)
    ## plt.plot(plArray, coeffs[2] + coeffs[1]*plArray+coeffs[0]*plArray**2)
    ## plt.title('FFT residuals by SPCTRANS position')
    ## plt.xlabel('SPCTRANS')
    ## plt.ylabel('Median of Residuals Between Top and Bottom Halves of FFT')

    # find the minimum; set the HPC path length position accordingly
    y_series = coeffs[2] + coeffs[1]*plArray+coeffs[0]*plArray**2
    max_y = min(y_series)  # find the minimum y-value
    zero_pl = pl_array[y_series.index(max_y)] # find the x-value (pathlength) corresponding to min y-value (residuals)
    # is the curve concave up?
    if (coeffs[0] < 0):
	print("Line fit to residuals is not concave up!!")

    # write data to file for copying and plotting on local machine
    file_name = "resids_test.csv"
    resids_df = pd.DataFrame(np.transpose([plArray,testArray]), columns=["pl","resid"])
    no_return = resids_df.to_csv(file_name)
    raw_input("Plot "+file_name+" on a local machine to see if OPD of zero has been found. Then proceed with [Enter]. Otherwise, break.")

    # command the HPC to move to that minimum
    pi.setINDI("Acromag.HPC.Tip=0;Tilt=0;Piston="+'{0:.1f}'.format(max_x)+";Mode=1")

    # as last step, remove grism
    raw_input("PRESS ENTER AFTER REMOVING GRISM AND INSERTING OBSERVING FILTERS")

    # turn off fizeau flag to avoid problems with other observations
    print("De-activating ROI aquisition flag")
    pi.setINDI("LMIRCAM.fizRun.value=Off")

    return


def optimize_opd_fizeau_airy(mode = "science"):
    # this dials OPD until the center of the coherence envelope is found

    # scan in OPD until there is a clear *global* maximum in the FFT_amp high-freq lobe amplitudes (i.e., the visibility of the fringes is highest)

    ## ## 1. take FFT of science PSF
    AmpPE, ArgPE = fft_img.fft(img_before_padding_before_FT)

    ## ## 2. measure amplitude of MTF high-freq node

    
    ## ## 3. command: move HPC in piston


    ## ## 4. repeat the above 3 steps, displaying the power each time and fitting a parabola to the data
    ## ## 5. let user press a key for either of the following two eventualities:
    ## ##    a. they are satisfied the fit is good, and command the HPC to go to the PL location of max power, or
    ## ##    b. the PL scan seems to be going in the wrong direction, and command the HPC to go back to the starting point and scan anew

    ## ## Last step: set FITS header flag 'FIZ_OPD_AIRY=1'
