from lmircam_tools import pi
import sys, os, string, time, pdb, copy, pyfits
import numpy as np
from pyindi import * 
from scipy import ndimage, sqrt, stats
#import matplotlib.pyplot as plt

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
# this is good as a fallback background subtraction method
# w = edge px on either side
def bkgdsub(image, method):

    ## ## found the below line; not sure what it was for
    ## ##image[image==0] = numpy.nan

    ctr = int( np.floor( image.shape[0] / 2. ) ) # half-width

    tmpimg = copy.copy(image) # make intermediary image

    # take mean or median across rows (for removing horizontal striping of array)
    if method == 'mean':
        #rowbkgd = stats.nanmean(tmpimg,1) # if deprecated, use 
        rowbkgd = np.nanmean(tmpimg,1)
    elif method == 'median':
        rowbkgd = np.nanmedian(tmpimg,1)
        rowbkgd2d = np.tile(np.reshape(rowbkgd,[len(rowbkgd),1]),[1,image.shape[1]]) # paint this column of median values into a 2048x2048 array
    tmpimg = np.subtract(tmpimg,rowbkgd2d) # simpl

    # do same as above, but for the columns
    if method == 'mean':
        colbkgd = np.nanmean(tmpimg,0)
    elif method == 'median':
        colbkgd = np.nanmedian(tmpimg,0)
        colbkgd2d = np.tile(np.reshape(colbkgd,[1,len(colbkgd)]),[image.shape[0],1])
        image = np.subtract(np.subtract(image,rowbkgd2d),colbkgd2d)

    image[np.isnan(image)] = 0

    # image now should have a nice flat background, but bad pixels will remain
    return image


# process the image (background-subtract and bad pix fix)
## ## ADD A FLAG TO ENABLE QUICK-AND-DIRTY BACKGROUND SUBTRACTION, WHICH HAS BEEN LEFT IN FOR NOW (EVEN THOUGH WE MAY BE SUBTRACTING ALREADY FROM THE WEB GUI)
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


