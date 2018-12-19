import sys, os, string, time, pdb, copy, pyfits
import numpy as np
#import np.ma as ma
from pyindi import *
import scipy
from scipy import ndimage, sqrt, stats, misc, signal
import pyfits
from lmircam_tools import *
from lmircam_tools import process_readout


def optimize_opd_fizeau_grism(psf_location, mode = "science"):
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

    # loop over some frames, and each time take FFT and move FPC in one direction
    for t in range(0,20):

	take_roi_background()

        print("Taking a background-subtracted frame")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=$

	if (mode == "fake_fits"):
            f = pyfits.open("test_frame_grismFiz_small.fits")

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

	# take right-hand side of FFT, integrate in x for detecting asymmetries
        to_corr = np.sum(AmpPE[:,int(0.5*img_before_padding_before_FT.shape[1]):], axis=1) 
        to_corr_masked = np.copy(to_corr)
	# mask the strong low-frequency power
	to_corr_masked[int(0.5*len(to_corr_masked))-2:int(0.5*len(to_corr_masked))+2] = 0 
	# take cross-correlation
        test_symm = signal.correlate(to_corr_masked, to_corr_masked[::-1], mode='same') 
	pdb.set_trace()

        leftHalf = test_symm[:int(0.5*len(test_symm))]
        rightHalf = test_symm[int(0.5*len(test_symm)):]
        resid = leftHalf-rightHalf[::-1] # find residuals of the cross-correlation sides
	pdb.set_trace()

        #plt.plot(resid)
        #frameArray = np.concatenate((frameArray,[f]))
        testArray = np.concatenate((testArray,[np.median(resid)])) # median of residuals array
        plArray = np.concatenate((plArray,[int(header['SPCTRANS'])])) # pathlength array
        #plt.show()
	pdb.set_trace()

        #plt.savefig("images/psf_altair_corr2_"+str("{:0>6d}".format(f))+".png", overwrite=False)
        #plt.clf()

	# now move the HPC
	stepSize = 5. # (um, total OPD)
	# big steps, translation stage: Ubcs.SPC_Trans.command=>5
	## ## pi.setINDI("Ubcs.SPC_Trans.command=>"+'{0:.1f}'.format(10*0.5*stepSize)) # factor of 10 bcz command is in 0.1 um
	# small steps, piezos: Acromag.HPC.Tip=0;Tilt=0;Piston=[stepSize];Mode=1
	## ## pi.setINDI("Acromag.HPC.Tip=0;Tilt=0;Piston="+'{0:.1f}'.format(stepSize)+";Mode=1")

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
    max_x = x[y_series.index(max_y)] # find the x-value corresponding to min y-value

    #########################################################

    # as last step, remove grism
    raw_input("PRESS ENTER AFTER REMOVING GRISM AND INSERTING OBSERVING FILTERS")


def optimize_opd_fizeau_airy(psf_location):
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
