import numpy as np
import matplotlib.pyplot as plt
import scipy
import numpy.ma as ma
from scipy import misc, signal
import astropy
from astropy.io import fits
from matplotlib.colors import LogNorm
import os.path

from lmircam_tools import *
from lmircam_tools import overlap_psfs

# AC Her
stem = ('/home/../../media/unasemaje/Seagate Expansion Drive/lbti_data_reduction/180524_fizeau_ac_her/raw_data/180524/')

# FOR 20180524 AC HER DATA
# for loop to make thumbnails of the data

#plt.figure()

# read in a simple background
str_backgrnd = stem+'lm_180524_'+str("{:0>6d}".format(2100))+'.fits'
backgd, header = fits.getdata(str_backgrnd,
                                 0,
                                 header=True) 

for f in range(9100,9102): #11497): # full dataset: (83,67168)

    filename_str = stem+'lm_180524_'+str("{:0>6d}".format(f))+'.fits'
    
    if os.path.isfile(filename_str): # if FITS file exists in the first place

        image, header = fits.getdata(filename_str,
                                 0,
                                 header=True) 
    
        # subt backgrnd and mask bullet bad parts
        #image = np.subtract(image,backgd)
        image[:,1100:1300] = np.median(image)
        image[:,:300] = np.median(image)
        image[:,-300:] = np.median(image)
        
        # locate PSF
        psf_loc = overlap_psfs.find_airy_psf(image)
    
        # make plot
        plt.imshow(image, origin="lower", cmap="gray")  
        #plt.ylabel('linear')
    
        # indicate PSF position
        plt.scatter([psf_loc[1],psf_loc[1]],[psf_loc[0],psf_loc[0]], 
                    color='red',s=90, edgecolors='r',facecolors='none') # put do
    
        plt.suptitle('(y,x) = '+str(psf_loc))
    
        #plt.show()
        
        # append PSF positions to file
        outF = open("ac_her_psf_coords.csv", "a")
        outF.write(str("{:0>6d}".format(f))+','+str(psf_loc[1])+','+str(psf_loc[0]))
        outF.write("\n")
    
        plt.tight_layout()
    
        plt.savefig("images/thumbnail_acher_180524_data_"+str("{:0>6d}".format(f))+".png", dpi=150, overwrite=False)
        plt.clf()
    
        #print('Frame '+str("{:0>6d}".format(f))+' done...')
outF.close()
