import sys, os, string, time, pdb, copy, pyfits
import numpy as np
import np.ma as ma
from pyindi import *
import scipy
from scipy import ndimage, sqrt, stats, misc, signal
import matplotlib.pyplot as plt
from astropy.io import fits
from lmircam_tools import * #process_readout


def optimize_tt_fizeau_airy(psf_location):
    # this changes the tip-tilt setpoints until the Fizeau Airy PSF looks good

    ## ## change FPC TT until phase gradients in PTF are removed
    ## ## 1. Differential tip: phase gradient is all up-down, and the low-freq node in FT amplitude takes on a crushed ellipticity.  
    ## ## 2. Differentia tilt: phase gradient is left-right, but it is not continuous– it is divided among the three nodes. 

    ## ## Last step: set FITS header flag 'FIZ_TT_AIRY=1'
