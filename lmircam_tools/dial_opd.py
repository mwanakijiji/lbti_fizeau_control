import sys, os, string, time, pdb, copy, pyfits
import numpy as np
import np.ma as ma
from pyindi import *
import scipy
from scipy import ndimage, sqrt, stats, misc, signal
import matplotlib.pyplot as plt
from astropy.io import fits
from lmircam_tools import * #process_readout

# this was first tested in testing_dial_opd_grism.ipynb

def optimize_opd_fizeau_grism(psf_location):

    ##'image' is from the science detector-- FILL IN!

    # cut out the grism image
    img_before_padding_before_FT = image[center_grism[0]-int(0.5*length_y):center_grism[0]+int(0.5*length_y),
                                     center_grism[1]-2*sig:center_grism[1]+2*sig]
    
    
    # take FFT
    
    # no padding for now
    AmpPE, ArgPE = fft_img.fft(img_before_padding_before_FT)

    #############################################################
    ## ## BEGIN OPTIMIZATION OF PATHLENGTH BASED ON GRISM FFTS
    plt.suptitle(str("{:0>6d}".format(f)))
    
    # testing
    to_corr = np.sum(AmpPE[:,int(0.5*img_before_padding_before_FT.shape[1]):], axis=1)
    to_corr_masked = np.copy(to_corr)
    to_corr_masked[98:102] = 0
    test_symm = signal.correlate(to_corr_masked, to_corr_masked[::-1], mode='same')
    
    leftHalf = test_symm[:int(0.5*len(test_symm))]
    rightHalf = test_symm[int(0.5*len(test_symm)):]
    resid = leftHalf-rightHalf[::-1]
    plt.plot(resid)
    frameArray = np.concatenate((frameArray,[f]))
    testArray = np.concatenate((testArray,[np.median(resid)]))
    plArray = np.concatenate((plArray,[int(header['SPCTRANS'])]))
    #plt.show()
    
    #plt.savefig("images/psf_altair_corr2_"+str("{:0>6d}".format(f))+".png", overwrite=False)
    plt.clf()

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
    ## ## END OPTIMIZATION OF PATHLENGTH
    #########################################################
    
    # as last step, remove grism


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
