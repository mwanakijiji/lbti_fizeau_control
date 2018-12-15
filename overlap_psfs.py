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


def find_grism_psf(image):
    if autoFindStar:

        # generate the Gaussian to correlate with the image
        mu_x_center = 0.5*image.shape[1] # initially center the probe shape for correlation
        sig = 5 # ersatz for now
	length_y = 40 # ersatz for now
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


def overlap_airy_psfs(psf_loc_setpoint):

    # take a background
    # allow aquisition from ROI box (keep smaller than 512x512!)
    ## ## PASSIVE MODE pi.setINDI("LMIRCAM.fizRun.value=On")
    print("Moving in a blank to take a background")
    ## ## PASSIVE MODE pi.setINDI("Lmir.lmir_FW4.command", "Blank", wait=True)
    print("Taking a background")
    ## ## PASSIVE MODE f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=1;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)


    ### move in half-moon to see SX first
    ## ## PASSIVE MODE pi.setINDI("Lmir.lmir_FW2.command", 'SX-Half-moon', timeout=45, wait=True)


    ### iterate to try to get SX PSF on the same pixel
    while True: #check_prereq_loops():

        ### locate SX PSF

	## ## NEED TO PUT IN A BLANK FOR A BACKGROUND

	# obtain new frame from LMIR
	print("Setting ROI aquisition flag")
	# allow aquisition from ROI box (keep smaller than 512x512!)
	## ## PASSIVE MODE pi.setINDI("LMIRCAM.fizRun.value=On")

	
	# take a frame with background subtracting
	print("Taking a background-subtracted frame")
	## ## PASSIVE MODE f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)

	## ## BEGIN TEST
	## ## read in fake FITS file
	#f = pyfits.open("test_frame.fits")
	## ## END TEST

	imgb4 = f[0].data
	imgb4bk = process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

	
	# locate the PSF
	## ## DO I REALLY WANT TO USE FIND_GRISM?
        psf_loc = find_grism_psf(imgb4bk)
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
	## ## PASSIVE MODE pi.setINDI("Acromag.FPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")

	
        ### re-locate SX PSF; correction needed?
        ## ## PASSIVE MODE pi.setINDI("LMIRCAM.fizRun.value=On")
        # take a frame with background subtracting
        print("Taking a background-subtracted frame")
        ## ## PASSIVE MODE f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)
	imgb4 = f[0].data
        imgb4bk = process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

	
	# locate the PSF
        ## ## DO I REALLY WANT TO USE FIND_GRISM?
        psf_loc = find_grism_psf(imgb4bk)

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
    ## ## PASSIVE MODE pi.setINDI("Lmir.lmir_FW2.command", 'DX-Half-moon', wait=True)


    ### NEW FOR-LOOP HERE
    while True: # do three iterations to try to get SX PSF on the same pixel

        # locate DX PSF
        ## ## PASSIVE MODE pi.setINDI("LMIRCAM.fizRun.value=On")
        # take a background
        print("Taking a background")
        ## ## PASSIVE MODE f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)
        # take a frame with background subtracting
        print("Taking a frame to analyze")
        ## ## PASSIVE MODE f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)

        imgb4bk =process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image
        psf_loc = find_airy_psf(imgb4bk) # locate the PSF

        ### MOVE HPC IN ONE STEP TO MOVE PSF TO RIGHT LOCATION
        vector_move_pix = np.subtract(psf_loc_setpoint,psf_loc) # vector of required movement in pixel space
        vector_move_asec = np.multiply(vector_move_pix,0.0107) # convert to asec
        ## ## PASSIVE MODE pi.setINDI("Acromag.HPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")

        ### RE-LOCATE DX PSF; CORRECTION NEEDED?
        ## ## PASSIVE MODE pi.setINDI("LMIRCAM.fizRun.value=On")
        # take a background
        print("Taking a background")
        ## ## PASSIVE MODE f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)
        # take a frame with background subtracting
        print("Taking a frame to analyze")
        ## ## PASSIVE MODE f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)

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
    ## ## PASSIVE MODE pi.setINDI("Lmir.lmir_FW2.command", 'Open', wait=True)




def overlap_grism_psfs(psf_loc_setpoint):
    # note psf_loc_setpoint is the desired location of the CENTER of the grism PSF

    # set some approximate parameters of the observed grism PSF
    sig = 5 # sigma of Gaussian profile in x (in pix)
    length_y = 200 # length in y of the psf (in pix)

    # define a Gaussian for correlating it to a grism-ed PSF in x dimension

    ### move in half-moon to see SX first
    ## ## PASSIVE MODE pi.setINDI("Lmir.lmir_FW2.command", 'SX-Half-moon', wait=True)#, timeout=45, wait=True)

    ### iterate to try to get SX PSF on the same pixel
    while True: 

        ### locate SX PSF
	print("Turning on fizRun")
	## ## PASSIVE MODE pi.setINDI("LMIRCAM.fizRun.value=On")
        ## ## PASSIVE MODE f=pi.getFITS("LMIRCAM.DisplayImage.File", "LMIRCAM.GetDisplayImage.Now", wait=True) # get what LMIR is seeing
        imgb4 = f[0].data
        #imgb4bk = bkgdsub(imgb4,'median') # simple background smoothing
        #imgb4bk -= numpy.median(imgb4bk) # subtract residual background

        imgb4bk =process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image

        psf_loc = find_grism_psf(imgb4bk) # locate the grism PSF center

        ### move FPC in one step to move PSF to right location
        vector_move_pix = np.subtract(psf_loc_setpoint,psf_loc) # vector of required movement in pixel space
        vector_move_asec = np.multiply(vector_move_pix,0.0107) # convert to asec
        ## ## PASSIVE MODE pi.setINDI("Acromag.FPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")

        ### re-locate SX PSF; correction needed?
        ## ## PASSIVE MODE f=pi.getFITS("LMIRCAM.DisplayImage.File", "LMIRCAM.GetDisplayImage.Now", wait=True) # get what LMIR is seeing
        imgb4 = f[0].data
        imgb4bk =process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image
        psf_loc = find_grism_psf(imgb4bk) # locate the PSF (y,x)
 
        print('-------------------')
        print('PSF location setpoint:')
        print(psf_loc_setpoint) 
        print('Current PSF loc:') 
        print(psf_loc)
    
        if (dist_pix(psf_loc,psf_loc_setpoint) < 5.):
            print('-------------------')
            print('Done moving one side. Switching to the other side.')
            break 

        print('Moving PSF again...')
    
    
    ### MOVE IN HALF-MOON TO SEE DX NEXT
    ## ## PASSIVE MODE pi.setINDI("Lmir.lmir_FW2.command", 'DX-Half-moon', timeout=45, wait=True)
    
    ### NEW FOR-LOOP HERE
    while True: # do three iterations to try to get SX PSF on the same pixel
    
        # locate DX PSF
        ## ## PASSIVE MODE f=pi.getFITS("LMIRCAM.DisplayImage.File", "LMIRCAM.GetDisplayImage.Now", wait=True) # get what LMIR is seeing
        imgb4 = f[0].data
                        
        imgb4bk =process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image
        psf_loc = find_grism_psf(imgb4bk) # locate the PSF
                                
        ### MOVE HPC IN ONE STEP TO MOVE PSF TO RIGHT LOCATION
        vector_move_pix = np.subtract(psf_loc_setpoint,psf_loc) # vector of required movement in pixel space
        vector_move_asec = np.multiply(vector_move_pix,0.0107) # convert to asec
        ## ## PASSIVE MODE pi.setINDI("Acromag.HPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")
                                                
        ### RE-LOCATE SX PSF; CORRECTION NEEDED?
        ## ## PASSIVE MODE f=pi.getFITS("LMIRCAM.DisplayImage.File", "LMIRCAM.GetDisplayImage.Now", wait=True) # get what LMIR is seeing
        imgb4 = f[0].data
        imgb4bk =process_readout.processImg(imgb4, 'median') # return background-subtracted, bad-pix-corrected image
        psf_loc = find_grism_psf(imgb4bk) # locate the PSF
    
        print('-------------------')
        print('PSF location setpoint:')
        print(psf_loc_setpoint) 
        print('Current DX PSF loc:')
        print(psf_loc)
                                                                            
        if (dist_pix(psf_loc,psf_loc_setpoint) < 5.):
            print('-------------------')
            print('Done moving DX. Switching to the other side.')
            break
                                                                                                    
        print('Moving DX again...')    


    print('Done moving PSFs. Reopening LMIR FW2.')
    ## ## PASSIVE MODE pi.setINDI("Lmir.lmir_FW2.command", 'Open', wait=True)
