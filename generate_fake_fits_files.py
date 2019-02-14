# This generates blank FITS files at a cadence meant to mimic LMIRcam during Fizeau observations.
# (These FITS files are read in by another script to extract FFT information.)

# created 2019 Jan 29 by E.S.


import pyfits
import numpy as np
import sched, time

retrieve_dir = "/opt/local/LBTI_INDI/data/LMIRCAM/180507/" # directory where we retrieve already-written frames
check_dir = "fake_monitor/" # directory for this night, which we are monitoring for new frames


def write_fake_fits(framenum):
    '''
    Write entirely new frames
    '''
    fake_frame = np.random.rand(2048,2048).astype(np.uint16)
    #hdu = pyfits.PrimaryHDU(fake_frame)
    #hdulist = pyfits.HDUList([hdu])
    hdu.writeto(check_dir+"lm_180507fake_"+str("{:0>6d}".format(framenum))+".fits", clobber=True)



def move_mimic_fits(framenum):
    '''
    Move old frames to another directory, to mimic
    writing new files to that directory
    '''
    hdulist = pyfits.open(retrieve_dir+"lm_180507_"+str("{:0>6d}".format(framenum))+".fits")
    hdulist.writeto(check_dir+"lm_180507_"+str("{:0>6d}".format(framenum))+".fits", clobber=True)



framenum = 8849

while True:
    time_start = time.time()
    time.sleep(1.0)
    #write_fake_fits(framenum)
    move_mimic_fits(framenum)
    framenum += 1
    print('written, time elapsed')
    print(str(time.time() - time_start))

