import numpy as np
import time
from numpy import ma
from pyindi import *
from scipy import ndimage, sqrt, stats, misc, signal

pi = PyINDI(verbose=False)

integ_time = 100 # integration time, msec

dir_to_monitor = "fake_monitor/" # directory to monitor for newly-written images
del_t = 1.0 # pause cadence (sec) with which to monitor that directory

#psf_loc_fake = (114,864) # (y,x) approximate location of PSF when reading in fake/old FITS files
psf_loc_fake = (292,853) # a grism sequence

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
        ArgPE = np.multiply(np.angle(PhaseExtract),180./np.pi)

        if mask:
            # mask out low-power regions
            AmpPE_masked = ma.masked_where(AmpPE < mask_thresh, AmpPE, copy=False)
            ArgPE_masked = ma.masked_where(AmpPE < mask_thresh, ArgPE, copy=False)
            return AmpPE_masked, ArgPE_masked

        else:
            return AmpPE, ArgPE


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


def take_roi_background(mode):
    ''' 
    This flips the Fizeau PSF flag to on, and takes a new background for that ROI
    '''

    # take a background
    # allow aquisition from ROI box (keep smaller than 512x512!)
    if (mode != "total_passive"):
        start_time = time.time()
        print("Setting ROI aquisition flag")
        pi.setINDI("LMIRCAM.fizRun.value=On")
        print("Moving in a blank to take a background")
        ## ## BELOW LINE WAS COMMENTED OUT DUE TO FW 4 STUCK 2018 dEC 29
        pi.setINDI("Lmir.lmir_FW4.command", "Blank", wait=True)
        print("Taking a background")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=1;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=60)
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


    if (image and (mode != "total_passive")):

        # turn on fizeau flag
        print("Activating ROI aquisition flag")
        pi.setINDI("LMIRCAM.fizRun.value=On")

        # take a new frame to see what things look like now
        pi.setINDI("LMIRCAM_save.enable_save.value=On")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=60)

        # turn off fizeau flag to avoid problems with other observations
        print("De-activating ROI aquisition flag")
        pi.setINDI("LMIRCAM.fizRun.value=Off")

        end_time = time.time()
        print("Grism put in in (secs)")
        print(end_time - start_time)
        print("-------------------")	

    return


def junk_empty_fcn(mode = "science"):
    print(mode)
    return
