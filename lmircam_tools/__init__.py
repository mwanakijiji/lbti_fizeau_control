import numpy as np
import pyfits
import time
from numpy import ma
from scipy import ndimage, sqrt, stats, misc, signal

# necessary for importing pyindi
import sys
sys.path.append('./eckharts_python_path')
from pyindi import *


#########################
### BEGIN USER INPUTS ###
#########################

# science wavelength
wavel_lambda = 3.87e-6 # filter central wavel (m)

# physical psf location, if code is being used to read in written FITS files
#psf_loc_fake = [-999,-999] # THIS TRIGGERS GAUSSIAN FIT TO DETERMINE PSF CENTER
#psf_loc_fake = [50,50] # (y,x) in FITS coordinates
#psf_loc_fake = [1277,569]
psf_loc_fake = [292,853] # a grism sequence
# desired OR current psf location, if code is using ROI/INDI machinery
fiz_lmir_sweet_spot = [50,50] # (y,x) in ROI coordinates

# directory to monitor for newly-written images
dir_to_monitor = "fake_monitor/"
#dir_to_monitor = "/home/observer/data/180507/"
#dir_to_monitor = "/opt/local/LBTI_INDI/data/LMIRCAM/190420/"
#dir_to_monitor = "/opt/local/LBTI_INDI/data/LMIRCAM/junk/"

# default background-subtraction options
# 1. indi_ROI: using indi ROIs, which requires a new frame to be taken with a blank in LMIR FW4
# 2. quick_dirt: just reading in written FITS files and doing a quick-and-dirty background subtraction
# N.b. some functions may request other type of background-subtraction
bkgd_choice = "quick_dirt"

#######################
### END USER INPUTS ###
####################### 


pi = PyINDI(verbose=True) # host needs to be default (lbti-web) to control stuff other than the fizeau driver
pi_fiz = PyINDI(host="localhost", verbose=True) # for accessing the fizeau driver

integ_time = 100 # integration time, msec (probably not needed; this is from when I was explicitly requesting $
del_t = 0.1 # pause cadence (sec) with which to monitor that directory
plateScale_LMIR = 0.0107 # in asec/pix

def needed_tt_setpt_corrxn(alpha,PS,Nx,Ny):
    ''' 
    Calculates the needed correction to the tip-tilt setpoints of the FPC,
    given a gradient vector of the PTF
    (see LBTI Controlled Fizeau research log, 2019 Feb 15)

    INPUTS:
    alpha: a vector [a,b] containing the gradient (in radians per Fourier pixel) in x and y
    PS: LMIR plate scale (asec per LMIR pixel)
    Nx: number of samples in x on LMIR (i.e., length in x of the array to be FFTed)
    Ny: same as Nx, but for y

    OUTPUT:
    corrxn_tilt_tip: a vector [delx,dely] of required tilt (x) and tip (y) changes to current setpoints
    '''

    delta = 1 # sampling spacing is 1 LMIR pixel
    beta_x = -alpha[0]*PS*Nx*delta/np.pi # (asec; the minus sign is because alpha and beta have opposite sign)
    beta_y = -alpha[1]*PS*Ny*delta/np.pi # (asec; the minus sign is because alpha and beta have opposite sign)

    # return the opposite (in units mas); this is the [tilt,tip] setpoint correction
    return -np.multiply([beta_x,beta_y],1000)


# define 2D gaussian for fitting PSFs
def gaussian_x(x, mu, sig):
    shape_gaussian = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    return shape_gaussian


def find_airy_psf(image):
    if True:

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

def find_grism_psf(image, sig, length_y):
    if True:

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

class fft_img:
    # take FFT of a 2D image

    def __init__(self, image):
        self.image = image

    def fft(self, padding=int(0), pad_mode='constant', mask_thresh=1e-10, mask=True):

        padI = np.pad(self.image, padding, pad_mode)
	# arguments: image, pad size, pad mode, threshold for masking, mask flag
        padI = np.fft.fftshift(padI)
        PhaseExtract = np.fft.fft2(padI)
        PhaseExtract = np.fft.fftshift(PhaseExtract)
        AmpPE = np.absolute(PhaseExtract)
        #ArgPE = np.multiply(np.angle(PhaseExtract),180./np.pi) # degrees
        ArgPE = np.angle(PhaseExtract) # radians

        print("Size of array being FFTed:")
        print(np.shape(PhaseExtract))

        if mask:
            # mask out low-power regions
            AmpPE_masked = ma.masked_where(AmpPE < mask_thresh, AmpPE, copy=False)
            ArgPE_masked = ma.masked_where(AmpPE < mask_thresh, ArgPE, copy=False)
            return AmpPE_masked, ArgPE_masked

        else:
            return AmpPE, ArgPE


class fft_img_ersatz:
    # this does nothing to an image except massage it in to a form that mimics Python-made FFTs
    # this allows me to take a FFTW image and let the downstream code accept it

    def __init__(self, image):
        self.image = image

    def fft(self, padding=int(0), pad_mode='constant', mask_thresh=1e-10, mask=True):

        padI = np.pad(self.image, padding, pad_mode)
        # arguments: image, pad size, pad mode, threshold for masking, mask flag
        AmpPE = padI
        #ArgPE = np.multiply(np.angle(PhaseExtract),180./np.pi) # degrees
        ArgPE = padI # radians

        print("Size of array being FFTed:")
        print(np.shape(PhaseExtract))

        if mask:
            # mask out low-power regions
            AmpPE_masked_ersatz = ma.masked_where(AmpPE < mask_thresh, AmpPE, copy=False)
            ArgPE_masked_ersatz = ma.masked_where(AmpPE < mask_thresh, ArgPE, copy=False)
            return AmpPE_masked_ersatz, ArgPE_masked_ersatz

        else:
            return AmpPE_ersatz, ArgPE_ersatz

# define 2D gaussian for fitting PSFs
def gaussian_x(x, mu, sig):
    shape_gaussian = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    return shape_gaussian


# find star and return coordinates [y,x]
def find_airy_psf(image):
    if True:

        #imageThis = numpy.copy(image)

        ''' 
        if (PSFside == 'left'):
            imageThis[:,1024:-1] = 0
        elif (PSFside == 'right'):
            imageThis[:,0:1024] = 0
        '''

        image[np.isnan(image)] = np.nanmedian(image) # if there are NaNs, replace them with the median image v$
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


def take_roi_background(mode, bkgd_mode):
    ''' 
    This flips the Fizeau PSF flag to on, and takes a new background for that ROI
    (probably obsolete, since user is suppsed to already take a background --2019 Mar 23, E.S.)
    '''

    # take a background
    # allow aquisition from ROI box (keep smaller than 512x512!)
    print("N.b. ROI background subtraction has been deactivated so as to do quick-and-dirty background subtraction with Vanessas algorithm")

    if ((mode != "total_passive") and (bkgd_mode != "quick_dirt")):
        start_time = time.time()
        print("Setting ROI aquisition flag")
        pi_fiz.setINDI("fizeau.enable_run.value=On")
        print("Moving in a blank to take a background")
        pi.setINDI("Lmir.lmir_FW4.command", "Blank", wait=True)
        print("Taking a background")
        f = pi_fiz.getFITS("fizeau.roi_image.file", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=1;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=60)
        print("Took a new background in (secs)")
        end_time = time.time()
        print(end_time - start_time)
        print("-------------------")

    return


def put_in_grism(mode = "science", image = "yes"):
    ''' 
    Inserts the LMIR grism
    '''

    start_time = time.time()

    fw25_selection = "Lspec2.8-4.0"
    fw3_selection = "Lgrism6AR"
    if (mode != "total_passive"):
        print("Entering grism mode: putting in "+fw25_selection+" and "+fw3_selection)
        pi.setINDI("Lmir.lmir_FW3.command", fw3_selection, timeout=45, wait=True)
        pi.setINDI("Lmir.lmir_FW25.command", fw25_selection, timeout=45, wait=True) # blocks some extraneous light

        end_time = time.time()
        print("Grism put in in (secs)")
        print(end_time - start_time)
        print("---------------------")

    return


def junk_empty_fcn(mode = "science"):
    print(mode)
    return
