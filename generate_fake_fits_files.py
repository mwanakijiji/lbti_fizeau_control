# This generates blank FITS files at a cadence meant to mimic LMIRcam during Fizeau observations.
# (These FITS files are read in by another script to extract FFT information.)

# created 2019 Jan 29 by E.S.


import pyfits
import numpy as np
import sched, time

datestring = "190419"
#retrieve_dir = "/opt/local/LBTI_INDI/data/LMIRCAM/190419/" # directory where we retrieve already-written frames
#retrieve_dir = "/opt/local/LBTI_INDI/data/LMIRCAM/180507/junk/"
retrieve_dir = "fake_retrieve/"
deposit_dir = "fake_monitor/" # directory for this night, which we are monitoring for new frames
#deposit_dir = "/opt/local/LBTI_INDI/data/LMIRCAM/junk/"


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
    hdulist = pyfits.open(retrieve_dir+"lm_"+datestring+"_"+str("{:0>6d}".format(framenum))+".fits")
    hdulist.writeto(deposit_dir+"lm_"+datestring+"_"+str("{:0>6d}".format(framenum))+"_junk.fits", clobber=True)



def move_fits_simple(framenum):
    ''' 
    Move a single fits frame to another directory, over and over, under new names
    '''

    #file_name_stem = "chrom_mono_avgwavel_5000_opd_00200_tip_0000_tilt_0000_transl_000_PS_10"
    #file_name_stem = "chrom_mono_avgwavel_5000_opd_00000_tip_0000_tilt_0010_transl_000_PS_10"
    #file_name_stem = "chrom_mono_avgwavel_5000_opd_00000_tip_0000_tilt_0090_transl_000_PS_10"
    #file_name_stem = "half_um_test"
    file_name_stem = "chrom_mono_avgwavel_5000_opd_00000_tip_0000_tilt_0014_transl_000_PS_10"
    hdulist = pyfits.open(retrieve_dir+file_name_stem+".fits")
    hdulist.writeto(deposit_dir+file_name_stem+"_"+str("{:0>6d}".format(framenum))+".fits", clobber=True)


### SEQUENCES OF FRAMES

# good Fizeau-grism sequences from 180507 Altair data:
# frames, rough PSF location
# 4249-4258, (y,x)=(302,852)
# 5706-5816, (y,x)=(292,853)
# 6949-7084, (y,x)=(328,818)

# one good Fizeau-Airy sequence from 180507 Altair data:
# 8849-9175, (y,x)=(115,864)

# Fizeau-Airy sequences around a fringe jump:
# 8259-8265, (y,x)=(115,864) (bright inner jailbar fringe slightly to the right)
# 8266-8284, (y,x)=Ibid. (looks symmetrical)
# 8285-8289, (y,x)=Ibid. (fringe jump; one dark jailbar down center)

start_framenum = 22800
stop_framenum = 22999
framenum = np.copy(start_framenum)

'''
while (framenum < stop_framenum):
    time_start = time.time()
    time.sleep(0.5)
    #write_fake_fits(framenum)
    move_mimic_fits(framenum)
    framenum += 1
    print('written, time elapsed')
    print(str(time.time() - time_start))
'''

### THE SAME FRAME, OVER AND OVER
start_framenum = 0
framenum = np.copy(start_framenum)
while True:
    time_start = time.time()
    time.sleep(0.5)
    #move_mimic_fits(framenum)
    move_fits_simple(framenum)
    framenum += 1
    print('written, time elapsed')
    print(str(time.time() - time_start))
