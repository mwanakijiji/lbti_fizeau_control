import sys, os, string, time, pdb, copy, pyfits
import numpy as np
import numpy.ma as ma
import scipy
from scipy import ndimage, sqrt, misc, stats,signal
from lmircam_tools import *
from lmircam_tools import pi, process_readout, gaussian_x, find_airy_psf, find_grism_psf


###### TO DO: THE OVERLAP FUNCTIONS HAVE A LOT OF SHARED FUNCTIONALITY; MAKE A CLASS STRUCTURE!
# SEE https://jeffknupp.com/blog/2014/06/18/improve-your-python-python-classes-and-object-oriented-programming/
######################################################################################


def dist_pix(current,goal):

    dist = np.sqrt(np.power((current[1]-goal[1]),2) + np.power((current[0]-goal[0]),2) )

    return dist


def centroid_and_move(side, tolerance = 5, mode = "science", psf_type = "airy"):
    ''' 
    Find the PSF and move telescope side or UBC mirror to place it at the setpoint

    INPUTS:
    side: "left" or "right" (there is no "both") ("Do. Or do not. There is no try." --Yoda)
    tolerance: how close (in pixels) do Airy PSF have to be positioned to the Fizeau position setpoint?
    mode: "fake_fits": read in fake FITS files (but continue sending LMIR and mirror commands)
	  "artif_source": use an artificial source (either laser or pinhole)
	  "science": on-sky
    psf_type: "airy" or "grism"

    PREREQS:
    1.   ROI Aquisition flag: LMIRCAM.fizRun.value=On
    2.	 A background frame must have been taken for the given ROI.
    '''

    if (side == "left"):
	half_moon_filter = "SX-Half-moon"
	x_side = "SX"
    elif (side == "right"):
	half_moon_filter = "DX-Half-moon"
	x_side = "DX"
    else:
	print("Side not defined!")
	return

    if (psf_type == "airy"):
	sig = 5
	length_y = 5
    elif (psf_type == "grism"):
	sig = 5
	length_y = 200

    ### move in half-moon to see the Airy PSF of one side
    print("Putting in "+half_moon_filter+" to see "+x_side)
    pi.setINDI("Lmir.lmir_FW2.command", half_moon_filter, timeout=45, wait=True)

    ### iterate to try to get Airy PSF on the same pixel
    while True:

        # take a frame with background subtracting
        print("Taking a background-subtracted frame")
	f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)

	if ((mode == "fake_fits") and (psf_type == "airy")):
            #f = pyfits.open("test_frame_fiz_large.fits")
            f = pyfits.open("test_fits_files/test_frame_fiz_small.fits")
	elif ((mode == "fake_fits") and (psf_type == "grism")):
	    f = pyfits.open("test_fits_files/test_frame_grismFiz_small.fits")

        imgb4 = f[0].data
        imgb4bk = process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

        # locate the PSF
        ## ## DO I REALLY WANT TO USE FIND_GRISM?
        psf_loc = find_grism_psf(imgb4bk, sig = sig, length_y = length_y)
        print('-------------')
        print(x_side + " PSF located at (ROI coords)")
        print(psf_loc)

        ### figure out required movement of PSF to right location
        vector_move_pix = np.subtract(psf_loc_setpoint,psf_loc) # vector of required movement in pixel space
        vector_move_asec = np.multiply(vector_move_pix,0.0107) # convert to asec
        print(x_side + " required vector movement in pix:")
        print(vector_move_pix)
        print(x_side + " required vector movement in asec:")
        print(vector_move_asec)

        # do a first movement with the telescope itself
        if (mode == "science"):
	    print("Moving " + x_side + " telescope")
            nod(XOff=vector_move_asec[1],YOff=vector_move_asec[0],side=side,Coords="DETXY",Type="REL")
            wait4AORunning() # let AO close

        ### re-locate PSF; correction needed?
        print("Taking a background-subtracted frame")
	f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)

	if (test_mode == "fake_fits"):
	    f = pyfits.open("test_fits_files/test_frame_fiz_small.fits")

	imgb4 = f[0].data
        imgb4bk = process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

        ## ## DO I REALLY WANT TO USE FIND_GRISM?
        psf_loc = find_grism_psf(imgb4bk, sig = sig, length_y = length_y) # find PSF

        print("-------------------")
        print("Fizeau PSF location setpoint:")
        print(psf_loc_setpoint)
        print("Current " + x_side + " PSF loc:")
        print(psf_loc)

        # if PSFs are closer than N pixels from each other, break
        ## ## TOLERANCE ON SKY SHOULD BE N=5 OR BETTER
        N = tolerance
        if (dist_pix(psf_loc,psf_loc_setpoint) < tolerance):
            print("-------------------")
            print("Done moving " + x_side)
            break

	else:
	    # fine-tune with FPC or HPC
	    if (side == "left"):
            	print("Moving SX PSF again, now with FPC movement")
            	pi.setINDI("Acromag.FPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")
	    elif (side == "right"):
		print("Moving DX PSF again, now with HPC movement")
		pi.setINDI("Acromag.HPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")

	if (mode == "fake_fits"):
	    # need to break, because otherwise the FPC/HPC mirrors won't converge
	    break


def overlap_psfs(psf_loc_setpoint, mode = "science", psf_type = "airy"):

    take_roi_background()

    raw_input("User: remove the Blank in FW4, then press return when done")

    centroid_and_move(side = "left", mode = mode, psf_type = psf_type)

    take_roi_background()

    centroid_and_move(side = "right", mode = mode, psf_type = psf_type)

    print('Done moving PSFs. Reopening LMIR FW2.')
    pi.setINDI("Lmir.lmir_FW2.command", 'Open', wait=True)

    # turn off fizeau flag to avoid problems with other observations
    print("De-activating ROI aquisition flag")
    pi.setINDI("LMIRCAM.fizRun.value=Off")

    return
