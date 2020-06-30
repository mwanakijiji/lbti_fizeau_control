import numpy as np
import pyfits
import time
from numpy import ma
from scipy import ndimage, sqrt, stats, misc, signal
from regions import PixCoord, CircleSkyRegion, CirclePixelRegion, PolygonPixelRegion


#################################################
# USER INPUTS: set units of the observation

D = 8.25 # diameter of one downstopped primary mirror (meters)
B_c2c = 14.4 # center-to-center separation (meters)
B_e2e = 22.65 # edge-to-edge separation (meters)
# science wavelength
wavel_lambda = 3.87e-6 # filter central wavel (m)
sci_wavel = wavel_lambda
plateScale_LMIR = 0.0107 # in asec/pix
plateScale = plateScale_LMIR # (asec/pixel)
global asecInRad # define global value
asecInRad = 206264.81 # (asec/rad)
#################################################

# physical psf location, if code is being used to read in written FITS files
#psf_loc_fake = [1024,1024] # (y,x) in FITS coordinates
#psf_loc_fake = [292,853] # a grism sequence
# desired OR current psf location, if code is using ROI/INDI machinery
#fiz_lmir_sweet_spot = [200,100] # (y,x) in ROI coordinates

# directory to monitor for newly-written images
dir_to_monitor = "fake_monitor/"
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
        # arguments: image, pad size, pad mode, threshold for masking, mask flag

        padI = np.pad(self.image, padding, pad_mode)
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

def fftMask(sciImg,wavel_lambda,plateScale,fyi_string=''):
    '''
    Take a FFT image, generate masks to select interesting areas of the FFT, and
    return data about those areas (amplitudes, normal vectors, etc.)

    INPUTS:
    sciImg: this is actually the FFT image, not the science detector image
    wavel_lambda: wavelength of the observation
    plateScale: plate scale of the detector (asec/pixel)
    fyi_string: an FYI string that could be used for plots

    OUTPUTS:
    dictFFTstuff: dictionary with keys corresponding to different parts of the FFT
    '''

    # make division lines separating different parts of the PSF
    line_M1diam_pixOnFFT = findFFTloc(D,np.shape(sciImg)[0],wavel_lambda,plateScale)
    line_center2center_pixOnFFT = findFFTloc(B_c2c,np.shape(sciImg)[0],wavel_lambda,plateScale)
    line_edge2edge_pixOnFFT = findFFTloc(B_e2e,np.shape(sciImg)[0],wavel_lambda,plateScale)

    # define circles
    circRad = 30 # pixels in FFT space
    circle_highFreqPerfect_L = CirclePixelRegion(center=PixCoord(x=line_center2center_pixOnFFT[0], y=0.5*np.shape(sciImg)[0]), radius=circRad)
    circle_highFreqPerfect_R = CirclePixelRegion(center=PixCoord(x=line_center2center_pixOnFFT[1], y=0.5*np.shape(sciImg)[0]), radius=circRad)
    circle_lowFreqPerfect = CirclePixelRegion(center=PixCoord(x=0.5*np.shape(sciImg)[1], y=0.5*np.shape(sciImg)[0]), radius=circRad)

    # define central rectangular region that includes all three nodes
    rect_pix = PolygonPixelRegion(vertices=PixCoord(x=[line_edge2edge_pixOnFFT[0],line_edge2edge_pixOnFFT[1],line_edge2edge_pixOnFFT[1],line_edge2edge_pixOnFFT[0]],
                                                       y=[line_M1diam_pixOnFFT[1],line_M1diam_pixOnFFT[1],line_M1diam_pixOnFFT[0],line_M1diam_pixOnFFT[0]]))

    # make the masks
    mask_circHighFreq_L = circle_highFreqPerfect_L.to_mask()
    mask_circHighFreq_R = circle_highFreqPerfect_R.to_mask()
    mask_circLowFreq = circle_lowFreqPerfect.to_mask()
    mask_rect = rect_pix.to_mask()

    ## apply the masks

    # initialize arrays of same size as science image
    sciImg1 = np.copy(sciImg)
    sciImg2 = np.copy(sciImg)
    sciImg3 = np.copy(sciImg)
    sciImg4 = np.copy(sciImg)
    sciImg5 = np.copy(sciImg)

    # region 1: high-freq lobe, left
    sciImg1.fill(np.nan) # initialize arrays of nans
    mask_circHighFreq_L.data[mask_circHighFreq_L.data == 0] = np.nan    # make zeros within mask cutout (but not in the mask itself) nans
    sciImg1[mask_circHighFreq_L.bbox.slices] = mask_circHighFreq_L.data  # place the mask cutout (consisting only of 1s) onto the array of nans
    sciImg1 = np.multiply(sciImg1,sciImg) # 'transmit' the original science image through the mask
    sciImg1 = sciImg1.filled(fill_value=np.nan) # turn all masked '--' elements to nans

    # region 2: high-freq lobe, right
    sciImg2.fill(np.nan) # initialize arrays of nans
    mask_circHighFreq_R.data[mask_circHighFreq_R.data == 0] = np.nan    # make zeros within mask cutout (but not in the mask itself) nans
    sciImg2[mask_circHighFreq_R.bbox.slices] = mask_circHighFreq_R.data  # place the mask cutout (consisting only of 1s) onto the array of nans
    sciImg2 = np.multiply(sciImg2,sciImg) # 'transmit' the original science image through the mask
    sciImg2 = sciImg2.filled(fill_value=np.nan) # turn all masked '--' elements to nans

    # region 3: low-freq lobe
    sciImg3.fill(np.nan) # initialize arrays of nans
    mask_circLowFreq.data[mask_circLowFreq.data == 0] = np.nan    # make zeros within mask cutout (but not in the mask itself) nans
    sciImg3[mask_circLowFreq.bbox.slices] = mask_circLowFreq.data  # place the mask cutout (consisting only of 1s) onto the array of nans
    sciImg3 = np.multiply(sciImg3,sciImg) # 'transmit' the original science image through the mask
    sciImg3 = sciImg3.filled(fill_value=np.nan) # turn all masked '--' elements to nans

    # region 4: rectangular region containing parts of all lobes
    sciImg4.fill(np.nan) # initialize arrays of nans
    mask_rect.data[mask_rect.data == 0] = np.nan    # make zeros within mask cutout (but not in the mask itself) nans
    sciImg4[mask_rect.bbox.slices] = mask_rect.data  # place the mask cutout (consisting only of 1s) onto the array of nans
    sciImg4 = np.multiply(sciImg4,sciImg) # 'transmit' the original science image through the mask
    sciImg4 = sciImg4.filled(fill_value=np.nan) # turn all masked '--' elements to nans
    import ipdb; ipdb.set_trace()

    # region 5: all three circular lobes
    sciImg5.fill(np.nan) # initialize arrays of nans
    sciImg5[mask_circLowFreq.bbox.slices] = mask_circLowFreq.data # low-freq lobe
    sciImg5[mask_circHighFreq_L.bbox.slices] = mask_circHighFreq_L.data # left high-freq lobe
    sciImg5[mask_circHighFreq_R.bbox.slices] = mask_circHighFreq_R.data # right high-freq lobe
    sciImg5 = np.multiply(sciImg5,sciImg) # 'transmit' the original science image through the mask
    sciImg5 = sciImg5.filled(fill_value=np.nan) # turn all masked '--' elements to nans
    import matplotlib.pyplot as plt
    plt.imshow(sciImg5,origin="lower")
    plt.savefig("junk_sciImg5.png")
    plt.close()
    plt.imshow(sciImg4,origin="lower")
    plt.savefig("junk_sciImg4.png")
    plt.close()

    # return medians of regions under masks
    med_highFreqPerfect_L = np.nanmedian(sciImg1)
    med_highFreqPerfect_R = np.nanmedian(sciImg2)
    med_lowFreqPerfect = np.nanmedian(sciImg3)
    med_rect = np.nanmedian(sciImg4)
    med_3_lobe = np.nanmedian(sciImg5)

    # return normal vectors corresponding to [x,y,z] to surfaces (x- and y- components are of interest)
    normVec_highFreqPerfect_L = normalVector(sciImg1)
    normVec_highFreqPerfect_R = normalVector(sciImg2)
    normVec_lowFreqPerfect = normalVector(sciImg3)
    normVec_rect = normalVector(sciImg4)
    normVec_3_lobe = normalVector(med_3_lobe)
    normVec_highFreqPerfect_L_x = normVec_highFreqPerfect_L[0]
    normVec_highFreqPerfect_L_y = normVec_highFreqPerfect_L[1]
    normVec_highFreqPerfect_L_z = normVec_highFreqPerfect_L[2]
    normVec_highFreqPerfect_R_x = normVec_highFreqPerfect_R[0]
    normVec_highFreqPerfect_R_y = normVec_highFreqPerfect_R[1]
    normVec_highFreqPerfect_R_z = normVec_highFreqPerfect_R[2]
    normVec_lowFreqPerfect_x = normVec_lowFreqPerfect[0]
    normVec_lowFreqPerfect_y = normVec_lowFreqPerfect[1]
    normVec_lowFreqPerfect_z = normVec_lowFreqPerfect[2]
    normVec_rect_x = normVec_rect[0]
    normVec_rect_y = normVec_rect[1]
    normVec_rect_z = normVec_rect[2]
    normVec_3_lobe_x = normVec_3_lobe[0]
    normVec_3_lobe_y = normVec_3_lobe[1]
    normVec_3_lobe_z = normVec_3_lobe[2]

    # return stdev in each region
    std_highFreqPerfect_L = np.nanstd(sciImg1)
    std_highFreqPerfect_R = np.nanstd(sciImg2)
    std_lowFreqPerfect = np.nanstd(sciImg3)
    std_rect = np.nanstd(sciImg4)
    std_3_lobe = np.nanstd(sciImg5)

    # generate images showing footprints of regions of interest
    # (comment this bit in/out as desired)

    # initialize dictionary to contain FFT data
    # N.b. all the info in this dictionary is EITHER for
    # the FFT amplitude OR the FFT phase, depending on what
    # the 'sciImg' is
    dictFFTstuff = {}

    # median of high-freq lobe on left side, within circular region centered around
    # where a perfect high-freq lobe would be
    dictFFTstuff["med_highFreqPerfect_L"] = med_highFreqPerfect_L

    # median of right-side high-freq lobe
    dictFFTstuff["med_highFreqPerfect_R"] = med_highFreqPerfect_R

    # median of low-frequency lobe
    dictFFTstuff["med_lowFreqPerfect"] = med_lowFreqPerfect

    # median of rectangle that is drawn to contain both high- and low-freq lobes
    dictFFTstuff["med_rect"] = med_rect

    # median of 3-lobe region containing both high- and low-freq lobes
    dictFFTstuff["med_three_lobe"] = med_3_lobe

    # stdev of the same regions
    dictFFTstuff["std_highFreqPerfect_L"] = std_highFreqPerfect_L

    # stdev of right-side high-freq lobe
    dictFFTstuff["std_highFreqPerfect_R"] = std_highFreqPerfect_R

    # stdev of low-frequency lobe
    dictFFTstuff["std_lowFreqPerfect"] = std_lowFreqPerfect

    # stdev of rectangle that is drawn to contain both high- and low-freq lobes
    dictFFTstuff["std_rect"] = std_rect

    # stdev of rectangle that is drawn to contain both high- and low-freq lobes
    dictFFTstuff["std_three_lobe"] = std_3_lobe

    # normal vectors to the high- and low- frequency
    # note vectors are [a,b,c] corresponding to the eqn Z = a*X + b*Y + c
    dictFFTstuff["dummy"] = [0.1] # dont touch! this lets me convert to a dataframe later
    dictFFTstuff["normVec_highFreqPerfect_L_x"] = normVec_highFreqPerfect_L_x
    dictFFTstuff["normVec_highFreqPerfect_L_y"] = normVec_highFreqPerfect_L_y
    dictFFTstuff["normVec_highFreqPerfect_L_z"] = normVec_highFreqPerfect_L_z
    dictFFTstuff["normVec_highFreqPerfect_R_x"] = normVec_highFreqPerfect_R_x
    dictFFTstuff["normVec_highFreqPerfect_R_y"] = normVec_highFreqPerfect_R_y
    dictFFTstuff["normVec_highFreqPerfect_R_z"] = normVec_highFreqPerfect_R_z
    dictFFTstuff["normVec_lowFreqPerfect_x"] = normVec_lowFreqPerfect_x
    dictFFTstuff["normVec_lowFreqPerfect_y"] = normVec_lowFreqPerfect_y
    dictFFTstuff["normVec_lowFreqPerfect_z"] = normVec_lowFreqPerfect_z
    dictFFTstuff["normVec_rect_x"] = normVec_rect_x
    dictFFTstuff["normVec_rect_y"] = normVec_rect_y
    dictFFTstuff["normVec_rect_z"] = normVec_rect_z
    dictFFTstuff["normVec_3_lobe_x"] = normVec_3_lobe_x
    dictFFTstuff["normVec_3_lobe_y"] = normVec_3_lobe_y
    dictFFTstuff["normVec_3_lobe_z"] = normVec_3_lobe_z

    # return the regions with mask overlays, too
    dictFFTstuff["sciImg1"] = sciImg1
    dictFFTstuff["sciImg2"] = sciImg2
    dictFFTstuff["sciImg3"] = sciImg3
    dictFFTstuff["sciImg4"] = sciImg4
    dictFFTstuff["sciImg5"] = sciImg5

    return dictFFTstuff


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


def findFFTloc(baseline,imageShapeAlong1Axis,wavel_lambda,plateScale,lOverD=1.):
    '''
    Returns the FFT pixel locations equivalent to a certain pixel distance on the science image
    '''

    # baseline: distance in physical space in the pupil plane (m)
    # imageShapeAlong1Axis: length of one side of the input image (needs to be square)
    # wavel_lambda: wavelength of light (m)
    # plateScale: detector plate scale (asec/pix)
    # lOverD: option if we are interested in the circular Airy rings (values 1.22, etc.); otherwise, just leave as 1

    line_diam_pixOnScience = lOverD*(wavel_lambda*asecInRad)/(baseline*plateScale) # distance in pixels on science detector
    line_diam_freq = np.divide(1.,line_diam_pixOnScience) # the corresponding frequency

    # find number of 'frequency pixels' multiply by units of l/D (like 1.22), and add to central x frequency pixel
    line_diam_pixOnFFT_L = 0.5*imageShapeAlong1Axis - np.divide(line_diam_freq,np.divide(1.,imageShapeAlong1Axis))
    line_diam_pixOnFFT_H = 0.5*imageShapeAlong1Axis + np.divide(line_diam_freq,np.divide(1.,imageShapeAlong1Axis))

    return line_diam_pixOnFFT_L, line_diam_pixOnFFT_H # the lower and higher values around freq of zero

def normalVector(sciImg):
    '''
    Fit a plane by finding the (normalized) normal vector to the best-fit plane

    INPUTS:
    sciImg: 2D image to fit a plane to

    OUTPUTS:
    normVec: the normal vector to the image data (note the normal vector is not normalIZED)
    '''

    # make x, y coords and flatten x, y, and z for putting into least-squares matrix equation
    X,Y = np.meshgrid(np.arange(0,np.shape(sciImg)[0]),
                  np.arange(0,np.shape(sciImg)[1]))
    Xflat_T = np.transpose(X.flatten())
    Yflat_T = np.transpose(Y.flatten())
    onesVec = np.ones(np.size(Xflat_T))
    Zflat_T = np.transpose(sciImg.flatten())

    # remove nans so we can fit planes
    Xflat_T = Xflat_T[~np.isnan(Zflat_T)]
    Yflat_T = Yflat_T[~np.isnan(Zflat_T)]
    onesVec = onesVec[~np.isnan(Zflat_T)]
    Zflat_T = Zflat_T[~np.isnan(Zflat_T)]

    # For a plane ax+by+c=z, the normal vector is [a,b,c]. To find this,
    # solve the matrix equation

    # AC = B, where
    # A are the x and y coords: [[x0,y0,1],[x1,y1,1],...[xN,yN,1]]
    # C are the coefficients we want: [[a],[b],[c]]
    # B is the z data: [[z0],[z1],...[zN]]

    # C = A+B, where A+ is the pseudoinverse of A, or A+ = ((A.T*A)^(-1))*A.T*B

    Amatrix = np.transpose(np.concatenate(([Xflat_T],[Yflat_T],[onesVec]), axis = 0))
    Bmatrix = Zflat_T

    # note that normVec = C
    normVec = np.dot(np.dot( np.linalg.pinv(np.dot(Amatrix.transpose(), Amatrix)), Amatrix.transpose()), Bmatrix)

    return normVec
