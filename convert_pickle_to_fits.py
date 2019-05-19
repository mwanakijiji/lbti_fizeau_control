# This notebook reads in simulated Fizeau PSFs and checks the gradients in the 
# PTF, so as to try to get a better intuition for what the gradient is saying 
# about needed tip-tilt offsets

# Created 2019 Feb 13 by E.S. (by branching off of (testing_making_set_masks.ipynb)

import numpy as np
import os.path
from astropy.io import fits
import pickle 
import glob

#from lmircam_tools import *
#from lmircam_tools import overlap_psfs

# make list of the available pickle files

file_list = glob.glob('./pickle_files/*pkl')


# if I want to read in a test frame of a simulated PSF, retrieve the pickle file and save the relevant
# array as a FITS file

# loop pickle files
for i in range(1,len(file_list)):

    print(file_list[i])

    '''
    name_stem = 'chrom_mono_avgwavel_5000_opd_00000_tip_'+\
                     str("{:0>4d}".format(int(tip_array[i])))+\
                     '_tilt_0000_transl_000_PS_10'
    '''
       
    dict_file = open(file_list[i], 'r')
    dict_pkld = pickle.load(dict_file)

    cube_init = np.zeros((3,np.shape(dict_pkld['PSF_image'])[0],np.shape(dict_pkld['PSF_image'])[1]))
    
    cube_init[0,:,:] = dict_pkld['PSF_image'].astype(np.float32)
    cube_init[1,:,:] = dict_pkld['FTamp_image'].astype(np.float32)
    cube_init[2,:,:] = dict_pkld['FTphase_image'].astype(np.float32)
    
    plt.imshow(dict_pkld['PSF_image'])
    plt.show()

    hdu = fits.PrimaryHDU(cube_init) # cubes with PSF, FFT_amp, FFT_arg
    #hdu = fits.PrimaryHDU(dict_pkld['PSF_image']) # single frame of PSF
    hdul = fits.HDUList([hdu])
    
    hdr = hdul[0].header
    hdr["OPD"] = dict_pkld['OPD']
    hdr["TIP_Y"] = dict_pkld['tip']
    hdr["TILT_X"] = dict_pkld['tilt']
    
    file_string = os.path.basename(file_list[i]).split(".")[0]
    
    hdul.writeto('./pickle_files/fits_equivalents/'+file_string+'.fits', overwrite=True)

