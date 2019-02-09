
# coding: utf-8

# In[ ]:

# This generates blank FITS files at a cadence meant to mimic LMIRcam during Fizeau observations.
# (These FITS files are read in by another script to extract FFT information.)

# created 2019 Jan 29 by E.S.


# In[1]:

import pyfits
import numpy as np
import sched, time


# In[2]:

retrieve_dir = "retrieve_dir" # directory where we retrieve already-written frames


# In[3]:

def write_fake_fits(framenum):
    '''
    Write entirely new frames
    '''
    fake_frame = np.random.rand(2048,2048).astype(np.uint16)
    hdu = pyfits.PrimaryHDU(fake_frame)
    #hdulist = pyfits.HDUList([hdu])
    hdu.writeto("check_dir/lm_180507fake_"+str("{:0>6d}".format(framenum))+".fits", clobber=True)


# In[5]:

def move_mimic_fits(framenum):
    '''
    Move old frames to another directory, to mimic
    writing new files to that directory
    '''
    fake_frame = pyfits.open("retrieve_dir/lm_180507_"+str("{:0>6d}".format(framenum))+".fits")
    hdu = pyfits.PrimaryHDU(fake_frame)
    #hdulist = pyfits.HDUList([hdu])
    hdu.writeto("check_dir/lm_180507fake_"+str("{:0>6d}".format(framenum))+".fits", clobber=True)


# In[4]:

framenum = 0

while True:
    time_start = time.time()
    time.sleep(0.2)
    #write_fake_fits(framenum)
    move_mimic_fits(framenum)
    framenum += 1
    print('written, time elapsed')
    print(str(time.time() - time_start))


# In[ ]:




# In[ ]:



