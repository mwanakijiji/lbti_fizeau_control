from lmircam_tools import pi
import sys, os, string, time, pdb, copy, pyfits
import numpy as np
from pyindi import * 
from scipy import ndimage, sqrt, stats
import matplotlib.pyplot as plt

# do simple for-loops to move the FPC and HPC

#pi.setINDI("LMIRCAM.Command.text=1 savedata")

#pi.setINDI("LMIRCAM.Command.text=go")
'''
for step in range(0,4):
    pi.setINDI("Acromag.HPC.Tip=0;Tilt=0.5;Piston=0;Mode=1")
    time.sleep(1.0)
    pi.setINDI("LMIRCAM.Command.text=go")
    pi.setINDI("Acromag.HPC.Tip=-0.5;Tilt=0.0;Piston=0;Mode=1")
    time.sleep(1.0)
    pi.setINDI("LMIRCAM.Command.text=go")
    pi.setINDI("Acromag.HPC.Tip=0;Tilt=-0.5;Piston=0;Mode=1")
    time.sleep(1.0)
    pi.setINDI("LMIRCAM.Command.text=go")
    pi.setINDI("Acromag.HPC.Tip=0.5;Tilt=0.0;Piston=0;Mode=1")
    time.sleep(1.0)
    pi.setINDI("LMIRCAM.Command.text=go")
    pi.setINDI("Acromag.FPC.Tip=0.0;Tilt=0.5;Piston=0;Mode=1")
    time.sleep(1.0)
    pi.setINDI("LMIRCAM.Command.text=go")
    pi.setINDI("Acromag.FPC.Tip=-0.5;Tilt=0;Piston=0;Mode=1")
    time.sleep(1.0)
    pi.setINDI("LMIRCAM.Command.text=go")
    pi.setINDI("Acromag.FPC.Tip=0.0;Tilt=-0.5;Piston=0;Mode=1")
    time.sleep(1.0)
    pi.setINDI("LMIRCAM.Command.text=go")
    pi.setINDI("Acromag.FPC.Tip=0.5;Tilt=0;Piston=0;Mode=1")
    time.sleep(1.0)
    pi.setINDI("LMIRCAM.Command.text=go")
'''

#####################################################
# START NEW CODE

# subtract mean or median (method) of the image
# w = edge px on either side
def bkgdsub(image, method):

    image[image==0] = numpy.nan

    ctr = int( np.floor( image.shape[0] / 2. ) ) # half-width

    tmpimg = copy.copy(image) # make intermediary image

    # take mean or median across rows (for removing horizontal striping of array)
    if method == 'mean':
        #rowbkgd = stats.nanmean(tmpimg,1) # if deprecated, use 
        rowbkgd = np.nanmean(tmpimg,1)
    elif method == 'median':
        rowbkgd = np.nanmedian(tmpimg,1)
    rowbkgd2d = np.tile(np.reshape(rowbkgd,[len(rowbkgd),1]),[1,image.shape[0]]) # paint this column of median values into a 2048x2048 array
    tmpimg = tmpimg - rowbkgd2d # simple form of background subtraction (vertical striping between columns will still be there)

    # do same as above, but for the columns
    if method == 'mean':
        colbkgd = np.nanmean(tmpimg,0)
    elif method == 'median':
        colbkgd = np.nanmedian(tmpimg,0)
    colbkgd2d = numpy.tile(np.reshape(colbkgd,[1,len(colbkgd)]),[image.shape[1],1])
    image = image - rowbkgd2d - colbkgd2d

    image[np.isnan(image)] = 0

    # image now should have a nice flat background, but bad pixels will remain
    return image


# process the image (background-subtract and bad pix fix)
def processImg(imgDummy, methodDummy):
    
    # bias level correction
    imgSub = bkgdsub(imgDummy,methodDummy) # simple background smoothing
    imgSub -= np.median(imgSub) # subtract residual background
    imgSubM = ndimage.median_filter(imgSub,3) # smoothed image

    # define BP mask
    imgDiff = numpy.abs(imgSub - imgSubM) # difference between smoothed and unsmoothed images
    stddev = numpy.std(imgDiff)
    mask = ( imgDiff > 4*stddev ) & ( imgDiff > 0.15 * numpy.abs(imgSub) ) # mask: True=bad, False=good

    imgSubBP = copy.copy(imgSub)
    imgSubBP[mask] = imgSubM[mask] # set bad pixels in unsmoothed image equal to those in smoothed image

    return imgSubBP


