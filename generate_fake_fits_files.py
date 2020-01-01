# This generates blank FITS files at a cadence meant to mimic LMIRcam during Fizeau observations.
# (These FITS files are read in by another script to extract FFT information.)

# created 2019 Jan 29 by E.S.


import pyfits
import numpy as np
import sched, time
import glob
import os

datestring = "180507"
#retrieve_dir = "/opt/local/LBTI_INDI/data/LMIRCAM/190419/" # directory where we retrieve already-written frames
#retrieve_dir = "/home/observer/data/180507/"
retrieve_dir = "fake_retrieve/synthetic/junk/"
#retrieve_dir = "fake_retrieve/"
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
    print(framenum)
    print(retrieve_dir+"lm_"+datestring)
    hdulist = pyfits.open(retrieve_dir+"lm_"+datestring+"_"+str("{:0>6d}".format(framenum))+".fits")
    hdulist.writeto(deposit_dir+"lm_"+datestring+"_"+str("{:0>6d}".format(framenum))+"_junk.fits", clobber=True)


def move_glob_fits(frame_string, csv_name):
    ''' 
    Move old frames to another directory, based on a string in their name
    '''
    fits_list = glob.glob(retrieve_dir + "*" + frame_string + "*.fits")
    fits_list_sorted = sorted(fits_list)

    print("File list is ")
    print(fits_list_sorted)

    # loop over each frame
    for i in range(0,len(fits_list_sorted)):
        hdulist = pyfits.open(fits_list_sorted[i])
        prihdr = hdulist[0].header
        hdulist.writeto(deposit_dir + os.path.basename(fits_list_sorted[i]), clobber=True)

        # meta-data
        opd = prihdr["OPD_UM"]
        tip = prihdr["TIPY_MAS"]
        tilt = prihdr["TILTXMAS"]

        # record the time when each frame is written, so I can compare true with retrieved values
        with open(csv_name, "a") as datalist:
            # cols of macie time; fake file name; OPD; tip; tilt
            datalist.write("%f, %s, %f, %f, %f\n" % (time.time(), os.path.basename(fits_list_sorted[i]), opd, tip, tilt))
            print("Deposited frame " + os.path.basename(fits_list_sorted[i]))
        time.sleep(0.14)



def move_fits_simple(framenum):
    ''' 
    Move a single fits frame to another directory, over and over, under new names
    '''
    framenum = int(framenum)
    #file_name_stem = "chrom_mono_avgwavel_5000_opd_00200_tip_0000_tilt_0000_transl_000_PS_10"
    #file_name_stem = "chrom_mono_avgwavel_5000_opd_00000_tip_0000_tilt_0010_transl_000_PS_10"
    #file_name_stem = "chrom_mono_avgwavel_5000_opd_00000_tip_0000_tilt_0090_transl_000_PS_10"
    #file_name_stem = "half_um_test"
    file_name_stem = "psf_trial4_00000921"
    hdulist = pyfits.open(retrieve_dir+file_name_stem+".fits")
    hdulist.writeto(deposit_dir+file_name_stem+"_"+str("{0:0>6d}".format(framenum))+".fits", clobber=True)


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


start_framenum = 5706
stop_framenum = 5816
framenum = np.copy(start_framenum)
'''
while (framenum < stop_framenum):
    time_start = time.time()
    time.sleep(1.0)
    #write_fake_fits(framenum)
    move_mimic_fits(framenum)
    framenum += 1
    print('written, time elapsed')
    print(str(time.time() - time_start))
'''

### THE SAME FRAME, OVER AND OVER
start_framenum = int(0)
framenum = np.copy(start_framenum)
while True:
    time_start = time.time()
    time.sleep(0.5)
    #move_mimic_fits(framenum)
    move_fits_simple(framenum)
    framenum += 1
    print('written, time elapsed')
    print(str(time.time() - time_start))

'''
### A GLOB OF FRAMES
frame_string = "trial1" # string for choosing the images
csv_string = frame_string + "_trial1_190612_injection.csv"
while True:
    time_start = time.time()
    move_glob_fits(frame_string, csv_string)
    print('written, time elapsed')
    print(str(time.time() - time_start))
'''
