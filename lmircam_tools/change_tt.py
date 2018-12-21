import numpy as np
import scipy
import numpy.ma as ma
import os.path
from scipy import misc, signal, ndimage
import pyfits
from astropy.coordinates import Angle, SkyCoord
from astropy.nddata.utils import extract_array
from regions import PixCoord, CircleSkyRegion, CirclePixelRegion, PolygonPixelRegion
from pyregion import read_region_as_imagecoord, get_mask
import csv
import time
import pickle
import pdb

from lmircam_tools import *
from lmircam_tools import overlap_psfs


#################################################
# USER INPUTS: set units of the observation

wavel_lambda = 4.051e-6 # filter central wavel (meters); n.b. Br-alpha is 4.051 um
D = 8.25 # diameter of one downstopped primary mirror (meters)
B_c2c = 14.4 # center-to-center separation (meters)
B_e2e = 22.65 # edge-to-edge separation (meters)
plateScale = 0.0107 # (asec/pixel)
global asecInRad # define global value
asecInRad = 206264.81 # (asec/rad)
#################################################


# locations of Airy minima/maxima (units lambda/D)
# (might use these later)
max1 = 0
min1 = 1.220
max2 = 1.635
min2 = 2.233
max3 = 2.679
min3 = 3.238
max4 = 3.699


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
    circRad = 60 # pixels in FFT space
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

    # return medians of regions under masks
    med_highFreqPerfect_L = np.nanmedian(sciImg1)
    med_highFreqPerfect_R = np.nanmedian(sciImg2)
    med_lowFreqPerfect = np.nanmedian(sciImg3)
    med_rect = np.nanmedian(sciImg4)

    # return normal vectors corresponding to [x,y,z] to surfaces (x- and y- components are of interest)
    normVec_highFreqPerfect_L = normalVector(sciImg1)
    normVec_highFreqPerfect_R = normalVector(sciImg2)
    normVec_lowFreqPerfect = normalVector(sciImg3)
    normVec_rect = normalVector(sciImg4)

    # return stdev in each region
    std_highFreqPerfect_L = np.nanstd(sciImg1)
    std_highFreqPerfect_R = np.nanstd(sciImg2)
    std_lowFreqPerfect = np.nanstd(sciImg3)
    std_rect = np.nanstd(sciImg4)

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

    # stdev of the same regions
    dictFFTstuff["std_highFreqPerfect_L"] = std_highFreqPerfect_L

    # median of right-side high-freq lobe
    dictFFTstuff["std_highFreqPerfect_R"] = std_highFreqPerfect_R

    # median of low-frequency lobe
    dictFFTstuff["std_lowFreqPerfect"] = std_lowFreqPerfect

    # median of rectangle that is drawn to contain both high- and low-freq lobes
    dictFFTstuff["std_rect"] = std_rect

    # normal vectors to the high- and low- frequency 
    # note vectors are [a,b,c] corresponding to the eqn Z = a*X + b*Y + c
    dictFFTstuff["normVec_highFreqPerfect_L"] = normVec_highFreqPerfect_L
    dictFFTstuff["normVec_highFreqPerfect_R"] = normVec_highFreqPerfect_R
    dictFFTstuff["normVec_lowFreqPerfect"] = normVec_lowFreqPerfect
    dictFFTstuff["normVec_rect"] = normVec_rect

    return dictFFTstuff


def get_pc_setpts(log_name = "setpt_log.csv", mode = "science"):
    ''' 
    Take FFT of PSF, and calculate new Phasecam PL and TT setpoints
    
    INPUTS:
    log_name: name of the csv file to which FFT information will be printed
    mode: "fake_fits": read in fake FITS files (but continue sending LMIR and mirror commands)
          "artif_source": use an artificial source (either laser or pinhole)
          "science": on-sky
    '''

    ampArray = []

    counter_num = 0

    take_roi_background()

    if (mode == "fake_fits"):
        #f = pyfits.open("test_fits_files/test_frame_fiz_large.fits")
        f = pyfits.open("test_fits_files/test_frame_fiz_small.fits")
    elif (mode == "science"):
        # take a frame with background subtracting
        print("Taking a background-subtracted frame")
        f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)

    image = f[0].data

    for f in range(0,1): # just 1 sample for now

        start = time.time() # start timer


        ## ##image, header = fits.getdata(filename_str,0,header=True)

        # locate PSF
        psf_loc = overlap_psfs.find_airy_psf(image)

        # size of cookie cut-out (measured center-to-edge)
        #cookie_size = 100 # maximum control radius as of 2018 July corresponds to 130.0 pixels

        # take FFT
        # cookie_cut = image[psf_loc[0]-cookie_size:psf_loc[0]+cookie_size,psf_loc[1]-cookie_size:psf_loc[1]+cookie_size]
        cookie_cut = np.copy(image)
        amp, arg = fft_img(cookie_cut).fft(padding=int(5*cookie_size), mask_thresh=1e5)

        # test: image with a perfect slope
        ''' 
        testing, header = fits.getdata('slope_test_psf.fits',0,header=True)
        cookie_cut_testing = testing[psf_loc[0]-cookie_size:psf_loc[0]+
                                     cookie_size,psf_loc[1]-cookie_size:psf_loc[1]+cookie_size]
        #sciImg = ma.asarray(sciImg)
        amp[np.isfinite(amp)] = -1 #cookie_cut_testing[np.isfinite(amp)]
        '''

        # sanity check (and to avoid getting for loop stuck)
        if (np.shape(amp)[0]!=np.shape(amp)[1]): # if the FFT doesn't make sense (i.e., if PSF was not found)
            print('PSF does not make sense ... aborting this one ...')
            continue

        print(amp.data)
        print(arg.data)

        # analyze FFTs
        fftInfo_amp = fftMask(amp,wavel_lambda,plateScale,
                                  fyi_string=str("{:0>6d}".format(f))+' FFT amp')
        fftInfo_arg = fftMask(arg,wavel_lambda,plateScale,
                                  fyi_string=str("{:0>6d}".format(f))+' FFT phase')

        print(fftInfo_amp)
        print(fftInfo_arg)

        # save fyi FITS files
        hdu = pyfits.PrimaryHDU(amp.data)
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto('junk_test_amp.fits', clobber=True)
        hdu = pyfits.PrimaryHDU(arg.data)
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto('junk_test_arg.fits', clobber=True)

        ## take info from the FFTs to send corrective movements

        # thresholds
        fft_ampl_high_freq_lowlimit = 2.4e5 # for good fringe visibility
        fft_ampl_low_freq_lowlimit = 1.4e6 # for acceptable AO correction
        fft_phase_vec_high_freq_highlimit = 5 # for Airy overlap
        std_lowFreqPerfect_lowlimit = 10 # for Airy overlap
        phase_normVec_highFreqPerfect_R_x # for phase of high-freq fringes

        # poor overlap of the Airy PSFs?
        print("--------------------------")
        print("Std of phase of low freq lobe:")
        print(fftInfo_arg["std_lowFreqPerfect"])
        ## HOW DOES IT COMPARE WITH std_lowFreqPerfect_lowlimit ?

        # poor fringe visibility
        print("Median of ampl of high freq lobe:")
        print(fftInfo_amp["med_highFreqPerfect_R"])
        ## HOW DOES IT COMPARE W fft_ampl_high_freq_lowlimit

        # high-freq fringes have strange phase
        ## ## change FPC TT until phase gradients in PTF are removed
        ## ## 1. Differential tip: phase gradient is all up-down, and the low-freq node in FT amplitude takes on a crushed ellipticity.
        ## ## 2. Differentia tilt: phase gradient is left-right, but it is not continuous; it is divided among the three nodes.
        print("Phase gradient in x of high freq in PTF:")
        print(fftInfo_arg["phase_normVec_highFreqPerfect_R"][0])
        print("Phase gradient in y of high freq in PTF:")
        print(fftInfo_arg["phase_normVec_highFreqPerfect_R"][1])

        # other quality control metrics from Phasecam (email from D. Defrere, 2018 Dec 17)
        # PCMSNR: S/N of K-band fringes
        # PCPHSTD: noise of phase in the integration time of NOMIC

        # all together now, lets make corrective movements
        # for better Airy overlap: tip-tilt the FPC
        pi.setINDI("Acromag.FPC.Tip="+'{0:.1f}'.format(vector_move_asec[0])+";Tilt="+'{0:.1f}'.format(vector_move_asec[1])+";Piston=0;Mode=1")

        # for better fringe visibility: move the FPC or HPC in piston
        stepSize = 5. # (um, total OPD)
        # big steps, translation stage: Ubcs.SPC_Trans.command=>5
        pi.setINDI("Ubcs.SPC_Trans.command=>"+'{0:.1f}'.format(10*0.5*stepSize)) # factor of 10 bcz command is in 0.1 um
        # small steps, piezos: Acromag.HPC.Tip=0;Tilt=0;Piston=[stepSize];Mode=1
        ## ## pi.setINDI("Acromag.HPC.Tip=0;Tilt=0;Piston="+'{0:.1f}'.format(stepSize)+";Mode=1")
        ## ## pi.setINDI("Acromag.FPC.Tip=0;Tilt=0;Piston="+'{0:.1f}'.format(stepSize)+";Mode=1")

        end = time.time()
        print(end - start)
        print('-----')

        # turn off fizeau flag to avoid problems with other observations
        print("De-activating ROI aquisition flag")
        pi.setINDI("LMIRCAM.fizRun.value=Off")

        return
