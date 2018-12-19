import sys, os, string, time, pdb, copy## ##, pyfits
import numpy as np
import numpy.ma as ma
## ##from pyindi import *
import scipy
from scipy import ndimage, sqrt, misc, stats,signal
#import matplotlib.pyplot as plt
import pyfits
from lmircam_tools import *
from lmircam_tools import process_readout
from lmircam_tools import pi

autoFindStar = True  # auto-detect star in frame?


###### TO DO: THE OVERLAP FUNCTIONS HAVE A LOT OF SHARED FUNCTIONALITY; MAKE A CLASS STRUCTURE!
# SEE https://jeffknupp.com/blog/2014/06/18/improve-your-python-python-classes-and-object-oriented-programming/
######################################################################################

# define 2D gaussian for fitting PSFs
## ## ONE COPY OF THIS FCN IS ALSO IN __INIT__; NEED TO CHANGE TO BE USING ONLY THE INIT VERSION
def gaussian_x(x, mu, sig):
    shape_gaussian = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    return shape_gaussian

# find star and return coordinates [y,x]
def find_airy_psf(image):
    if autoFindStar:

        #imageThis = numpy.copy(image)

        '''
        if (PSFside == 'left'):
            imageThis[:,1024:-1] = 0
        elif (PSFside == 'right'):
            imageThis[:,0:1024] = 0
        '''

        image[np.isnan(image)] = np.nanmedian(image) # if there are NaNs, replace them with the median image value
        imageG = ndimage.gaussian_filter(image, 6) # further remove effect of bad pixels (somewhat redundant?)
        loc = np.argwhere(imageG==imageG.max())
        cx = loc[0,1]
        cy = loc[0,0]

        #plt.imshow(imageG, origin="lower")
        #
        #plt.scatter([cx,cx],[cy,cy], c='r', s=50)
        #plt.colorbar()
        #plt.show()
        #print [cy, cx] # check

    return [cy, cx]


## ## THIS HAS ANOTHER COPY IN __INIT__; CHANGE TO ONLY USE __INIT__ COPY
def find_grism_psf(image, sig, length_y):
    if autoFindStar:

        # generate the Gaussian to correlate with the image
        mu_x_center = 0.5*image.shape[1] # initially center the probe shape for correlation
	x_probe_gaussian = gaussian_x(np.arange(image.shape[1]), mu_x_center, sig)

        # generate a top hat for correlating it to a grism-ed PSF in y
        y_abcissa = np.arange(image.shape[0])
        y_probe_tophat = np.zeros(image.shape[0])
        y_probe_tophat[np.logical_and(y_abcissa > 0.5*image.shape[0]-0.5*length_y,
                          y_abcissa <= 0.5*image.shape[0]+0.5*length_y)] = 1

        # dampen edge effects
        repl_val = np.median(image)
        image[0:4,:] = repl_val
        image[-5:-1,:] = repl_val
        image[:,0:4] = repl_val
        image[:,-5:-1] = repl_val

        # correlate real PSF and probe shape in x
        corr_x = signal.correlate(np.sum(image,axis=0), x_probe_gaussian, mode='same')

        # correlate real PSF and probe shape in y
        # (first, we have to integrate along x and subtract the background shape)
        profile_y = np.subtract(np.sum(image,axis=1), image.shape[1]*np.median(image,axis=1)) 
        corr_y = signal.correlate(profile_y, y_probe_tophat, mode='same')

        # find centers of psfs
        psf_center_x = np.argmax(corr_x)
        psf_center_y = np.argmax(corr_y)

    return [psf_center_y, psf_center_x]



def dist_pix(current,goal):

    dist = np.sqrt(np.power((current[1]-goal[1]),2) + np.power((current[0]-goal[0]),2) )

    return dist


def centroid_and_move(side, mode = "science"):
    ''' 
    Find the PSF and move telescope side or UBC mirror to place it at the setpoint

    INPUTS:
    side: "left" or "right" (there is no 'both')
    mode: "fake_fits": read in fake FITS files (but continue sending LMIR and mirror commands)
	  "artif_source": use an artificial source (either laser or pinhole)
	  "science": on-sky

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
	break

    ### move in half-moon to see the Airy PSF in question
    pi.setINDI("Lmir.lmir_FW2.command", half_moon_filter, timeout=45, wait=True)

    ### iterate to try to get Airy PSF on the same pixel
    while True:

        # take a frame with background subtracting
        print("Taking a background-subtracted frame")

	f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)

	if (mode == "fake_fits"):
            #f = pyfits.open("test_frame_fiz_large.fits")
            f = pyfits.open("test_frame_fiz_small.fits")

        imgb4 = f[0].data
        imgb4bk = process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

        # locate the PSF
        ## ## DO I REALLY WANT TO USE FIND_GRISM?
        psf_loc = find_grism_psf(imgb4bk, 5, 5)
        print('-------------')
        print(x_side + " PSF located at")
        print(psf_loc)

        ### figure out required movement of PSF to right location
        vector_move_pix = np.subtract(psf_loc_setpoint,psf_loc) # vector of required movement in pixel space
        vector_move_asec = np.multiply(vector_move_pix,0.0107) # convert to asec
        print(x_side + " vector movement in pix:")
        print(vector_move_pix)
        print(x_side + " vector movement in asec:")
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
	    f = pyfits.open("test_frame_fiz_small.fits")

	imgb4 = f[0].data
        imgb4bk = process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

        # locate the PSF
        ## ## DO I REALLY WANT TO USE FIND_GRISM?
        psf_loc = find_grism_psf(imgb4bk, 5, 5)

        print("-------------------")
        print("Fizeau PSF location setpoint:")
        print(psf_loc_setpoint)
        print("Current " + x_side + " PSF loc:")
        print(psf_loc)

        # if PSFs are closer than N pixels from each other, break
        ## ## TOLERANCE ON SKY SHOULD BE N=5
        N = 50.
        if (dist_pix(psf_loc,psf_loc_setpoint) < N):
            print("-------------------")
            print("Done moving " + x_side)
            break

	else:
	    # fine-tune with FPC or HPC
	    if (side == "left"):
            	print("Moving SX PSF again, now with FPC movement")
            	pi.setINDI("Acromag.FPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")
	    elif (side == "right"):
		pi.setINDI("Acromag.HPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")

	if (mode == "fake_fits"):
	    # need to break, because otherwise the FPC/HPC mirrors won't converge
	    break

    


def overlap_airy_psfs(psf_loc_setpoint):

    # take a background
    # allow aquisition from ROI box (keep smaller than 512x512!)
    print("Setting ROI aquisition flag")
    pi.setINDI("LMIRCAM.fizRun.value=On")
    print("Moving in a blank to take a background")
    pi.setINDI("Lmir.lmir_FW4.command", "Blank", wait=True)
    print("Taking a background")
    f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=1;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)

    raw_input("User: remove the Blank in FW4, then press return when done")

    centroid_and_move(side = "left")

    # take a background
    # allow aquisition from ROI box (keep smaller than 512x512!)
    pi.setINDI("LMIRCAM.fizRun.value=On")
    print("Moving in a blank to take a background")
    pi.setINDI("Lmir.lmir_FW4.command", "Blank", wait=True)
    print("Taking a background")
    f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=1;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)
    ## ## TAKE BLANK BACK OUT

    centroid_and_move(side = "right")

    print('Done moving PSFs. Reopening LMIR FW2.')
    pi.setINDI("Lmir.lmir_FW2.command", 'Open', wait=True)




def overlap_grism_psfs(psf_loc_setpoint):
    # note psf_loc_setpoint is the desired location of the CENTER of the grism PSF

    # set some approximate parameters of the observed grism PSF
    sig = 5 # sigma of Gaussian profile in x (in pix)
    length_y = 200 # length in y of the psf (in pix)

    # take a background
    # allow aquisition from ROI box (keep smaller than 512x512!)
    pi.setINDI("LMIRCAM.fizRun.value=On")
    print("Moving in a blank to take a background")
    pi.setINDI("Lmir.lmir_FW4.command", "Blank", wait=True)
    print("Taking a background")
    f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=1;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)


    ### move in half-moon to see SX first
    pi.setINDI("Lmir.lmir_FW2.command", 'SX-Half-moon', timeout=45, wait=True)


    ### iterate to try to get SX PSF on the same pixel
    while True: #check_prereq_loops():

        ### locate SX PSF

        # obtain new frame from LMIR
        print("Setting ROI aquisition flag")
        # allow aquisition from ROI box (keep smaller than 512x512!)
        pi.setINDI("LMIRCAM.fizRun.value=On")

        # take a frame with background subtracting
        print("Taking a background-subtracted frame")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=$

        ## ## BEGIN TEST
        ## ## read in fake FITS file
        #f = pyfits.open("test_frame_fiz_large.fits")
        f = pyfits.open("test_frame_grismFiz_small.fits")
        ## ## END TEST

        imgb4 = f[0].data
        imgb4bk = process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

        # locate the PSF
        ## ## DO I REALLY WANT TO USE FIND_GRISM?
        psf_loc = find_grism_psf(imgb4bk, sig, length_y)
        print('-------------')
        print('SX PSF located at')
        print(psf_loc)

        ### move FPC in one step to move PSF to right location
        vector_move_pix = np.subtract(psf_loc_setpoint,psf_loc) # vector of required movement in pixel space
        vector_move_asec = np.multiply(vector_move_pix,0.0107) # convert to asec
        print("SX vector movement in pix:")
        print(vector_move_pix)
        print("SX vector movement in asec:")
        print(vector_move_asec)
        print("Making FPC movement")
        pi.setINDI("Acromag.FPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")

        ### re-locate SX PSF; correction needed?
        pi.setINDI("LMIRCAM.fizRun.value=On")
        # take a frame with background subtracting
        print("Taking a background-subtracted frame")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=$
        imgb4 = f[0].data
        imgb4bk = process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

        # locate the PSF
        ## ## DO I REALLY WANT TO USE FIND_GRISM?
        psf_loc = find_grism_psf(imgb4bk, sig, length_y)

        print('-------------------')
        print('Fizeau PSF location setpoint:')
        print(psf_loc_setpoint)
        print('Current SX PSF loc:')
        print(psf_loc)

        # if PSFs are closer than N pixels from each other, break
        ## ## TOLERANCE ON SKY SHOULD BE N=5
        N = 50.
        if (dist_pix(psf_loc,psf_loc_setpoint) < N):
            print('-------------------')
            print('Done moving SX. Switching to DX.')
            break

        print('Moving SX PSF again...')


    ### MOVE IN HALF-MOON TO SEE DX NEXT
    pi.setINDI("Lmir.lmir_FW2.command", 'DX-Half-moon', wait=True)


    ### NEW FOR-LOOP HERE
    while True: # do three iterations to try to get SX PSF on the same pixel

        # locate DX PSF
        pi.setINDI("LMIRCAM.fizRun.value=On")
        # take a background
        print("Taking a background")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=$
        # take a frame with background subtracting
        print("Taking a frame to analyze")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=$

        imgb4bk =process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image
        psf_loc = find_airy_psf(imgb4bk) # locate the PSF

        ### MOVE HPC IN ONE STEP TO MOVE PSF TO RIGHT LOCATION
        vector_move_pix = np.subtract(psf_loc_setpoint,psf_loc) # vector of required movement in pixel space
        vector_move_asec = np.multiply(vector_move_pix,0.0107) # convert to asec
        pi.setINDI("Acromag.HPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")

        ### RE-LOCATE DX PSF; CORRECTION NEEDED?
        pi.setINDI("LMIRCAM.fizRun.value=On")
        # take a background
        print("Taking a background")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=$
        # take a frame with background subtracting
        print("Taking a frame to analyze")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=$

        print('-------------------')
        print('Fizeau PSF location setpoint:')
        print(psf_loc_setpoint)
        print('Current DX PSF loc:')
        print(psf_loc)

        ## ## TOLERANCE ON SKY SHOULD BE N=5
        N = 50
        if (dist_pix(psf_loc,psf_loc_setpoint) < 50.):
            print('-------------------')
            print('Done moving DX.')
            break

        print('Moving DX again...')


    print('Done moving PSFs. Reopening LMIR FW2.')
    pi.setINDI("Lmir.lmir_FW2.command", 'Open', wait=True)
