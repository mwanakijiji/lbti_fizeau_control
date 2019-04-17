import sys, os, string, time, pdb, copy, pyfits
import numpy as np
import numpy.ma as ma
import scipy
import time
from scipy import ndimage, sqrt, misc, stats,signal
from lmircam_tools import *
from lmircam_tools import pi, process_readout, gaussian_x, find_airy_psf, find_grism_psf
from lmircam_tools.utils import nod, wait4AORunning


###### TO DO: THE OVERLAP FUNCTIONS HAVE A LOT OF SHARED FUNCTIONALITY; MAKE A CLASS STRUCTURE!
# SEE https://jeffknupp.com/blog/2014/06/18/improve-your-python-python-classes-and-object-oriented-programming/
######################################################################################


def dist_pix(current,goal):

    dist = np.sqrt(np.power((current[1]-goal[1]),2) + np.power((current[0]-goal[0]),2) )

    return dist


def centroid_and_move(psf_loc_setpt, side, tolerance = 5, mode = "science", psf_type = "airy"):
    ''' 
    Find the PSF and move telescope side or UBC mirror to place it at the setpoint

    INPUTS:
    side: "left" or "right" (there is no "both") ("Do. Or do not. There is no try." --Yoda)
    tolerance: how close (in pixels) do Airy PSF have to be positioned to the Fizeau position setpoint?
    mode: testing or science
    psf_type: "airy" or "grism"

    PREREQS:
    ROI Aquisition flag: fizeau.enable_run.value=On
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

        # take a background-subtracted frame
        print("Taking a background-subtracted frame")
        pi.setINDI("LMIRCAM_save.enable_save.value=On")
        f = pi_fiz.getFITS("fizeau.roi_image.file", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=60)

        # get the right image slice
        if (np.ndim(f[0].data) > 2):
            imgb4 = f[0].data[-1,:,:] # images from LMIRcam (> summer 2018) are cubes of nondestructive reads
        else:
            imgb4 = np.squeeze(f[0].data)

        imgb4bk = np.copy(imgb4) # vestigial; image should already be background-subtracted

        # locate the PSF
        ## ## (do I really want to use find_grism for an Airy PSF? it seems to work pretty well...)
        psf_loc = find_grism_psf(imgb4bk, sig = sig, length_y = length_y)
        print('-------------')
        print(x_side + " PSF located at ROI coords (y,x) in pix")
        print(psf_loc)

        ### figure out required movement of PSF to right location
        vector_move_pix = np.subtract(psf_loc_setpt,psf_loc) # vector of required movement in pixel space
        vector_move_asec = np.multiply(vector_move_pix,0.0107) # convert to asec
        print(x_side + " required vector movement in (y,x) pix:")
        print(vector_move_pix)
        print(x_side + " required vector movement in (y,x) asec:")
        print(vector_move_asec)

        # do a first movement with the telescope itself
        if (mode == "science"):
	    print("Moving " + x_side + " telescope")
            nod(XOff=vector_move_asec[1],YOff=vector_move_asec[0],side=side,Coords="DETXY",Type="REL")
            wait4AORunning('both') # let AO close

        ### re-locate PSF; correction needed?
        print("Taking a background-subtracted frame")
        pi.setINDI("LMIRCAM_save.enable_save.value=On")
	f = pi_fiz.getFITS("fizeau.roi_image.file", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=60)

        # make sure we get right image slice
        if (np.ndim(f[0].data) > 2):
            imgb4 = f[0].data[-1,:,:] # images from LMIRcam (> summer 2018) are cubes of nondestructive reads
        else:
            imgb4 = np.squeeze(f[0].data)

        imgb4bk = np.copy(imgb4) # vestigial; image should already be background-subtracted

        psf_loc = find_grism_psf(imgb4bk, sig = sig, length_y = length_y) # find PSF

        print("-------------------")
        print("Fizeau PSF location setpoint in (y,x) pix:")
        print(psf_loc_setpt)
        print("Current " + x_side + " PSF loc in (y,x) pix:")
        print(psf_loc)

        # if PSFs are closer than N pixels from each other, break
        ## ## TOLERANCE ON SKY SHOULD BE N=5 OR SMALLER
        N = tolerance
        if (dist_pix(psf_loc,psf_loc_setpt) < tolerance):
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


def overlap_psfs(integ_time, fiz_lmir_sweet_spot, mode = "science", psf_type = "airy"):

    start_time = time.time()

    raw_input("User: remove the Blank in FW4, then press return when done")

    centroid_and_move(fiz_lmir_sweet_spot, side = "left", mode = mode, psf_type = psf_type)

    centroid_and_move(fiz_lmir_sweet_spot, side = "right", mode = mode, psf_type = psf_type)

    print('Done moving PSFs. Reopening LMIR FW2.')
    pi.setINDI("Lmir.lmir_FW2.command", 'Open', wait=True)

    # take a new frame to see what things look like now
    pi.setINDI("LMIRCAM_save.enable_save.value=On")
    f = pi_fiz.getFITS("fizeau.roi_image.file", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=60)

    # turn off fizeau flag to avoid problems with other observations
    print("De-activating ROI aquisition flag")
    pi_fiz.setINDI("fizeau.enable_run.value=Off")
    end_time = time.time()
    print("PSF overlapping done in (secs)")
    print(end_time - start_time)
    print("-------------------")

    return
