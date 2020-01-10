#!/usr/bin/env python
# coding: utf-8

# This notebook makes thumbnails of lots of Fizeau frames so as to isolate what
# I want to do engineering tests

# created from a parent, 2020 Jan. 9 by E.S.

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os.path
import glob

from lmircam_tools import process_readout

# get all filenames in the directory
dir_name = ('./')
file_names_all = glob.glob(dir_name + "*.fits")

for f in range(0,len(file_names_all)):
    
    image, header = fits.getdata(file_names_all[f],0,header=True) 
    
    if (np.ndim(image) > 2):
        image = image[-1,:,:] 
    else:
        image = np.squeeze(image)

    # a quick and dirty background subtraction and bad pixel fixing
    image = process_readout.processImg(image,"median",background=True) # simple background subtraction
    
    # make plot
    plt.imshow(image, origin="lower", cmap="gray") # PSF cut-out
    plt.ylabel('y')
    plt.xlabel('x')
    plt.suptitle(str(os.path.basename(file_names_all[f])))

    plt.tight_layout()
    plt.savefig("thumbnail_"+str(os.path.basename(file_names_all[f]).split(".")[0])+".png", dpi=100, overwrite=True)
    plt.clf()
    
    print('Frame '+str(os.path.basename(file_names_all[f]))+' done...')
    

