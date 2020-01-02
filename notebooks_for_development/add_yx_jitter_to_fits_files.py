# This reads in FITS files and translates them in (y,x) random walks, to simulate
# jitter and centroiding error

# Created 2020 Jan 1 by E.S. (from parent .ipynb notebook)

import numpy as np
#import matplotlib.pyplot as plt
import scipy
from scipy import ndimage
from astropy.io import fits
from walker_def import walker

x_scale = 25 # in pix
y_scale = 25 # in pix
N = 10000

xArray = walker(x_scale, N)
yArray = walker(y_scale, N)


# trial5 synthetic dataset: just take one frame and move it around

original_filename_str = "psf_trial1_00000001.fits"
for f in range(0,N):
    # read in same image
    image, hdr = fits.getdata(original_filename_str,0,header=True)

    x_shift = xArray[f]
    y_shift = yArray[f]

    # translate it
    image_translated = scipy.ndimage.interpolation.shift(image,[y_shift,x_shift])

    # update header
    hdr["X_SHIFT"] = x_shift
    hdr["Y_SHIFT"] = y_shift

    # write out
    write_file_name = "psf_trial5_" + str("{:0>8d}".format(int(f))) + ".fits"
    fits.writeto(write_file_name,
                 data=image_translated, header=hdr, overwrite=True)
    print("Wrote out " + str(write_file_name))


# trial6 synthetic dataset: take all frames, one at a time, from trial1 data
# so as to have opd, tip, tilt, and (y,x) random walks all at once

for f in range(0,N):

    # read in a new image each time
    original_filename_str = "psf_trial1_" + str("{:0>8d}".format(int(f))) + ".fits"
    image, hdr = fits.getdata(original_filename_str,0,header=True)

    x_shift = xArray[f]
    y_shift = yArray[f]

    # translate it
    image_translated = scipy.ndimage.interpolation.shift(image,[y_shift,x_shift])

    # update header
    hdr["X_SHIFT"] = x_shift
    hdr["Y_SHIFT"] = y_shift

    # write out
    write_file_name = "psf_trial6_" + str("{:0>8d}".format(int(f))) + ".fits"
    fits.writeto(write_file_name,
                 data=image_translated, header=hdr, overwrite=True)
    print("Wrote out " + str(write_file_name))
