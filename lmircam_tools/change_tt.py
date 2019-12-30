import numpy as np
import scipy
import numpy.ma as ma
import os.path
from scipy import misc, signal, ndimage
import pyfits
#from astropy.coordinates import Angle, SkyCoord
#from astropy.nddata.utils import extract_array
#from regions import PixCoord, CircleSkyRegion, CirclePixelRegion, PolygonPixelRegion
#from pyregion import read_region_as_imagecoord, get_mask
import csv
import time
import pickle
import pdb
import pandas as pd
import math
import glob

from lmircam_tools import *
from lmircam_tools import overlap_psfs
from lmircam_tools import process_readout

#################################################
# USER INPUTS: set units of the observation

D = 8.25 # diameter of one downstopped primary mirror (meters)
B_c2c = 14.4 # center-to-center separation (meters)
B_e2e = 22.65 # edge-to-edge separation (meters)
plateScale = plateScale_LMIR # (asec/pixel)
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
    circRad = 30 # pixels in FFT space

    '''
    ## Commented out because the machinery is not compatible with Python 2.6.6 on lbti-macie --E.S., 2019 Dec. 29
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
    '''
    # initialize arrays of same size as science image
    sciImg1 = np.copy(sciImg)
    sciImg2 = np.copy(sciImg)
    sciImg3 = np.copy(sciImg)
    sciImg4 = np.copy(sciImg)

    # region 1: high-freq lobe, left
    sciImg1.fill(np.nan) # initialize arrays of nans
    x1_center = line_center2center_pixOnFFT[0]
    y1_center = 0.5*np.shape(sciImg)[0]
    sciImg1[y1_center-circRad:y1_center+circRad, x1_center-circRad:x1_center+circRad] = 1.
    ##mask_circHighFreq_L.data[mask_circHighFreq_L.data == 0] = np.nan    # make zeros within mask cutout (but not in the mask itself) nans
    ##sciImg1[mask_circHighFreq_L.bbox.slices] = mask_circHighFreq_L.data  # place the mask cutout (consisting only of 1s) onto the array of nans
    sciImg1 = np.multiply(sciImg1,sciImg) # 'transmit' the original science image through the mask
    ##sciImg1 = sciImg1.filled(fill_value=np.nan) # turn all masked '--' elements to nans

    # region 2: high-freq lobe, right
    sciImg2.fill(np.nan) # initialize arrays of nans
    x2_center = line_center2center_pixOnFFT[1]
    y2_center = 0.5*np.shape(sciImg)[0]
    sciImg2[y2_center-circRad:y2_center+circRad, x2_center-circRad:x2_center+circRad] = 1.
    ##mask_circHighFreq_R.data[mask_circHighFreq_R.data == 0] = np.nan    # make zeros within mask cutout (but not in the mask itself) nans
    ##sciImg2[mask_circHighFreq_R.bbox.slices] = mask_circHighFreq_R.data  # place the mask cutout (consisting only of 1s) onto the array of nans
    sciImg2 = np.multiply(sciImg2,sciImg) # 'transmit' the original science image through the mask
    ##sciImg2 = sciImg2.filled(fill_value=np.nan) # turn all masked '--' elements to nans

    # region 3: low-freq lobe
    sciImg3.fill(np.nan) # initialize arrays of nans
    x3_center = 0.5*np.shape(sciImg)[1]
    y3_center = 0.5*np.shape(sciImg)[0]
    sciImg3[y3_center-circRad:y3_center+circRad, x3_center-circRad:x3_center+circRad] = 1.
    ##mask_circLowFreq.data[mask_circLowFreq.data == 0] = np.nan    # make zeros within mask cutout (but not in the mask itself) nans
    ##sciImg3[mask_circLowFreq.bbox.slices] = mask_circLowFreq.data  # place the mask cutout (consisting only of 1s) onto the array of nans
    sciImg3 = np.multiply(sciImg3,sciImg) # 'transmit' the original science image through the mask
    ##sciImg3 = sciImg3.filled(fill_value=np.nan) # turn all masked '--' elements to nans

    # region 4: rectangular region containing parts of all lobes
    sciImg4.fill(np.nan) # initialize arrays of nans
    sciImg4[line_M1diam_pixOnFFT[0]:line_M1diam_pixOnFFT[1], line_edge2edge_pixOnFFT[0]:line_edge2edge_pixOnFFT[1]] = 1.
    ##mask_rect.data[mask_rect.data == 0] = np.nan    # make zeros within mask cutout (but not in the mask itself) nans
    ##sciImg4[mask_rect.bbox.slices] = mask_rect.data  # place the mask cutout (consisting only of 1s) onto the array of nans
    sciImg4 = np.multiply(sciImg4,sciImg) # 'transmit' the original science image through the mask
    ##sciImg4 = sciImg4.filled(fill_value=np.nan) # turn all masked '--' elements to nans

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

    # stdev of right-side high-freq lobe
    dictFFTstuff["std_highFreqPerfect_R"] = std_highFreqPerfect_R

    # stdev of low-frequency lobe
    dictFFTstuff["std_lowFreqPerfect"] = std_lowFreqPerfect

    # stdev of rectangle that is drawn to contain both high- and low-freq lobes
    dictFFTstuff["std_rect"] = std_rect

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

    # return the regions with mask overlays, too
    dictFFTstuff["sciImg1"] = sciImg1
    dictFFTstuff["sciImg2"] = sciImg2
    dictFFTstuff["sciImg3"] = sciImg3
    dictFFTstuff["sciImg4"] = sciImg4

    return dictFFTstuff


def print_write_fft_info(integ_time, sci_wavel, mode = "science", setpoints_pickle = " ", checker=False):
    '''
    Take FFT of PSF, and calculate new Phasecam PL and TT setpoints

    INPUTS:
    mode: testing or science
    setpoints_pickle: filename of pickle file containing previous setpoints, so we can compare
    checker=True: this function is being run to check effect of a previous correction

    OUTPUTS:
    N/A; FFT info from series of Fizeau PSFs is pickled
    '''

    ampArray = []
    fftInfo_amp_df = pd.DataFrame()
    fftInfo_arg_df = pd.DataFrame()

    counter_num = 0 # for counting number of analyzed PSFs

    # read in any new images written out to a directory
    files_start = glob.glob(dir_to_monitor + "*.fits") # starting list of files
    num_psfs_to_analyze = 10 # number of PSFs to sample (a very large number if just writing retrieved values from single frames to file)

    while counter_num < num_psfs_to_analyze:

        # if we're just reading in fake FITS files, monitor a directory
        # and read in a new file
        if (mode == "fake_fits"):

            while True:

            	time_start = time.time()
            	time.sleep(del_t)

            	fft_pickle_write_name = "pickled_info/fft_info_"+str("{:0>2d}".format(counter_num))+".pkl" # filename for pickled FFT info
            	#if (checker == True): # if we are checking a previous correction, lets distinguish the pickle file name
            	#    fft_pickle_write_name = "fft_info_check_"+str("{:0>2d}".format(counter_num))+".pkl"

            	# check to see if there were new files from last check
            	files_later = glob.glob(dir_to_monitor + "/*.fits")

            	# are there new files?
            	new_list = np.setdiff1d(files_later,files_start)

            	time_start = time.time() # start timer

            	'''
            	# this is for individual test images
            	if ((mode == "fake_fits") or (mode == "total_passive")):
                    files_start = glob.glob(dir_to_monitor + "*.fits")
                    #f = pyfits.open("test_fits_files/test_frame_fiz_large.fits")
                    f = pyfits.open("test_fits_files/test_frame_fiz_small.fits")
                elif ((mode == "science") or (mode == "nac_source") or (mode == "az_source")):
                    # take a frame with background subtracting
                    print("Taking a background-subtracted frame")
                    pi.setINDI("lmircam_save.enable_save.value=On")
                f = pi_fiz.getFITS("fizeau.roi_image.file", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=60)
                '''

            	# if there are no new files, cycle back through
            	if (len(new_list) <= 2):
              	    continue

            	# if there are new files, read one of them in
            	elif (len(new_list) > 2):

                    print("New files detected.")

                    # reassign these files to be next starting point
                    files_start = files_later

                    # filename of second-newst file (in case the newest is still being written)
                    second_newest = sorted(files_later)[-2]

                    # read in file
                    f = pyfits.open(second_newest)
                    file_name_full = os.path.basename(second_newest) # filename without the path
                    file_name_base = os.path.splitext(file_name_full)[0] # filename without the extension, too

                    # break out of the 'while True'
                    break

    	# if we're reading in real frames from LMIR, catch them as they are read out
    	elif ((mode == "az_source") or (mode == "science")):

            time_start = time.time()

            file_name_full = str(int(time.time())) + "_" + str(counter_num) # string for written file names

            # roi image
            f = pi_fiz.getFITS("fizeau.roi_image.file", timeout=60)
            file_name_base = "test_placeholder" # string to add to log image files

            # FFT mag image (this has to wait for a new frame)
            fftw_amp = pi_fiz.getFITS("fizeau.mag_image.file", timeout=60)

            # FFT phase image (this has to wait for a new frame(
            fftw_phase = pi_fiz.getFITS("fizeau.phase_image.file", timeout=60)

    	# get the right image slice
        # if the frames are live images from LMIR, take the slice with the longest integration time
        if ((mode == "az_source") or (mode == "science")):
    	    if (np.ndim(f[0].data) > 2):
                image = f[0].data[-1,:,:] # images from LMIRcam (> summer 2018) are cubes of nondestructive reads
    	    else:
                image = np.squeeze(f[0].data)
        # if the frames are fake FITS files, I want the first slice
        elif (mode == "fake_fits"):
            if (np.ndim(f[0].data) > 2):
                image = f[0].data[0,:,:]
            else:
                image = np.squeeze(f[0].data)

    	# if this is a fake fits file, do a quick-and-dirty background subtraction
        if (mode == "fake_fits"):
            image = process_readout.processImg(image,"median") # simple background subtraction

    	# save detector image to check (overwrites previous)
    	hdu = pyfits.PrimaryHDU(image)
    	hdulist = pyfits.HDUList([hdu])
    	hdu.writeto("junk_test_image_seen.fits", clobber=True)

    	# locate PSF
    	psf_loc = find_grism_psf(image, sig=5, length_y=5)
    	if (mode == "fake_fits"):
            psf_loc = psf_loc_fake # if we are reading in fake FITS files
            ## ## if time, test the following line
            #psf_loc = find_airy_psf(image)
            print("Found PSF location:")
            print(psf_loc)

        # size of cookie cut-out (measured center-to-edge)
        cookie_size = 50 # maximum control radius as of 2018 July corresponds to 130.0 pixels

        # take FFT (some of this is redundant if INDI is putting out FFTs)
        if (mode == "fake_fits"):
            cookie_cut = image[psf_loc[0]-cookie_size:psf_loc[0]+cookie_size,psf_loc[1]-cookie_size:psf_loc[1]+cookie_size]
        else:
            cookie_cut = np.copy(image)
    	padding_choice = int(5*cookie_size)
        print("Total length of one side of image being FFTed (pix):")
        print(np.shape(cookie_cut)[0]+2*padding_choice)
        amp, arg = fft_img(cookie_cut).fft(padding=padding_choice, mask_thresh=1e5)

        # this is a kludge for slipping in the INDI FFT amplitude in place of the Python one
        # (the phase has a checkerboard pattern until Paul fixes it, so Im just going to keep
        # the Python phase)
        if ((mode == "az_source") or (mode == "science")):
            amp = fftw_amp[0]
            amp_ersatz, junk = fft_img(fftw_amp[0].data).fft(padding=padding_choice, mask_thresh=1e5)
            #amp = np.pad(amp, padding_choice, 'constant') # padding necessary to be consistent with Python FFT phase array size
            #amp = np.ma.asarray(amp) # turning it into a masked array necessary to feed it into fftMask()
            amp = amp_ersatz

    	# save image to check
    	hdu = pyfits.PrimaryHDU(cookie_cut)
    	hdulist = pyfits.HDUList([hdu])
    	hdu.writeto("log_images/cookie_prepadding_prefft_" + file_name_base + ".fits", clobber=True)
    	# test: see what the FFT looks like
        hdu = pyfits.PrimaryHDU(amp.data)
    	hdulist = pyfits.HDUList([hdu])
    	hdu.writeto("log_images/fft_amp_" + file_name_base + ".fits", clobber=True)
        hdu = pyfits.PrimaryHDU(arg.data)
    	hdulist = pyfits.HDUList([hdu])
    	hdu.writeto("log_images/fft_arg_" + file_name_base + ".fits", clobber=True)
    	# --commented out because it was triggering on NxM frames where N!=M--
    	# sanity check (and to avoid getting for loop stuck)
    	#if (np.shape(amp)[0] != np.shape(amp)[1]): # if the FFT doesn't make sense (i.e., if PSF was not found)
    	#    print('PSF does not make sense ... aborting this one ...')
    	#    continue

    	# analyze FFTs
    	fftInfo_amp = fftMask(amp,sci_wavel,plateScale,
                                  fyi_string=" FFT amp")
        fftInfo_arg = fftMask(arg,sci_wavel,plateScale,
                                  fyi_string=" FFT phase")

        ## write out MTF info for analysis later
        # take strip 5 pixels tall in the middle
        y_side_length = np.shape(amp.data)[0]
        strip = amp.data[int(0.5*y_side_length)-2:int(0.5*y_side_length)+3,:]
        med_strip = np.median(strip, axis = 0)
        max_strip = np.max(strip, axis = 0)
        # integrate the whole MTF in y (so that we can do a fuller integration later,
        # but still have data columns of the right length to write)
        y_integ = np.sum(amp.data, axis = 0)
        # put data into dictionary
        d = {"x_pix": np.arange(np.shape(amp.data)[1]),
             "med_strip": med_strip,
             "max_strip":max_strip,
             "y_integ":y_integ}
        # convert to df
        mtf_data = pd.DataFrame(data = d)
        # write to file
        mtf_data.to_csv("log_images/mtf_data/mtf_data_" + file_name_base + ".csv")

        print(type(np.ma.asarray(amp)))
        print(type(arg))
        print(type(amp.data))
        print(type(arg.data))
        print(np.shape(amp.data))
        print(np.shape(arg.data))
    	# save fyi FITS files to see the masks, etc.
    	hdu = pyfits.PrimaryHDU(fftInfo_amp["sciImg1"])
    	hdulist = pyfits.HDUList([hdu])
    	hdu.writeto("log_images/fft_amp_masked_region_1_" + file_name_base + ".fits", clobber=True)
    	hdu = pyfits.PrimaryHDU(fftInfo_amp["sciImg2"])
    	hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/fft_amp_masked_region_2_" + file_name_base + ".fits", clobber=True)
    	hdu = pyfits.PrimaryHDU(fftInfo_amp["sciImg3"])
    	hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/fft_amp_masked_region_3_" + file_name_base + ".fits", clobber=True)
    	hdu = pyfits.PrimaryHDU(fftInfo_amp["sciImg4"])
    	hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/fft_amp_masked_region_4_" + file_name_base + ".fits", clobber=True)
        hdu = pyfits.PrimaryHDU(fftInfo_arg["sciImg1"])
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/fft_arg_masked_region_1_" + file_name_base + ".fits", clobber=True)
        hdu = pyfits.PrimaryHDU(fftInfo_arg["sciImg2"])
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/fft_arg_masked_region_2_" + file_name_base + ".fits", clobber=True)
        hdu = pyfits.PrimaryHDU(fftInfo_arg["sciImg3"])
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/fft_arg_masked_region_3_" + file_name_base + ".fits", clobber=True)
        hdu = pyfits.PrimaryHDU(fftInfo_arg["sciImg4"])
        hdulist = pyfits.HDUList([hdu])
        hdu.writeto("log_images/fft_arg_masked_region_4_" + file_name_base + ".fits", clobber=True)
    	#hdu = pyfits.PrimaryHDU(amp.data)
    	#hdulist = pyfits.HDUList([hdu])
    	#hdu.writeto("junk_other_tests/junk_test_amp_psf.fits", clobber=True)
    	#hdu = pyfits.PrimaryHDU(arg.data)
    	#hdulist = pyfits.HDUList([hdu])
    	#hdu.writeto("junk_other_tests/junk_test_arg_psf.fits", clobber=True)

    	## take info from the FFTs to send corrective movements

    	# thresholds
    	fft_ampl_high_freq_lowlimit = 2.4e5 # min threshold for good fringe visibility
    	fft_ampl_low_freq_lowlimit = 1.4e6 # min threshold for acceptable AO correction
    	fft_phase_vec_high_freq_highlimit = 5 # max threshold for Airy overlap
    	std_lowFreqPerfect_lowlimit = 10 # max threshold for Airy overlap
    	phase_normVec_highFreqPerfect_R_x = 100 # max threshold for phase of high-freq fringes


    	## MONITOR REALTIME STATUS OF PSF QUALITY CRITERIA

    	# Overlap of the Airy PSFs?
        print("--------------------------")
    	print("--------------------------")
    	print("Std of phase of low freq lobe:")
    	print(fftInfo_arg["std_lowFreqPerfect"])
    	# TO CORRECT: TT THE HPC TO TRANSLATE DX AIRY PSF
    	## ## (MAYBE TRY TTING THE FPC LATER?)

    	# High-freq fringe visibility (median)
    	print("--------------------------")
    	print("Median of ampl of high freq lobe:")
    	print(fftInfo_amp["med_highFreqPerfect_R"])
    	# TO CORRECT: MOVE THE SPC_TRANS TO FIND CENTER OF COHERENCE ENVELOPE

    	# High-freq phase gradient
    	print("--------------------------")
        print("Median of phase (science PSF):")
        print(fftInfo_arg["med_highFreqPerfect_R"]*(180./np.pi))
        print("Pathlength correction needed (um):")
        print(fftInfo_arg["med_highFreqPerfect_R"]*((180./np.pi)/360.)*sci_wavel*1e6)
    	print("Phase gradient in x of high freq in PTF:")
    	print(fftInfo_arg["normVec_highFreqPerfect_R_x"])
    	print("Phase gradient in y of high freq in PTF:")
    	print(fftInfo_arg["normVec_highFreqPerfect_R_y"])
    	print("--------------------------")
    	# TO CORRECT: TT THE FPC? (NOT SURE HERE...)
    	## ## N.B. 1. Differential tip: phase gradient is all up-down, and the low-freq node in FT amplitude takes on a crushed ellipticity.
    	## ## N.B. 2. Differential tilt: phase gradient is left-right, but it is not continuous; it is divided among the three nodes.

    	## HOW DO THE ABOVE PRINTED QUANTITIES COMPARE WITH THE THRESHOLDS?

    	# other quality control metrics from Phasecam (email from D. Defrere, 2018 Dec 17)
    	# PCMSNR: S/N of K-band fringes
    	# PCPHSTD: noise of phase in the integration time of NOMIC

    	# also find time
    	timestamp = time.time()
    	fftInfo_amp["time"] = timestamp
    	fftInfo_arg["time"] = timestamp

    	# convert dictionaries to dataframes that are easy to write to file
    	# have to delete the 2D arrays from the dictionaries so they can be converted to dataframes
    	del fftInfo_amp["sciImg1"], fftInfo_amp["sciImg2"], fftInfo_amp["sciImg3"], fftInfo_amp["sciImg4"]
    	del fftInfo_arg["sciImg1"], fftInfo_arg["sciImg2"], fftInfo_arg["sciImg3"], fftInfo_arg["sciImg4"]
    	#fftInfo_amp_df = fftInfo_amp.drop(['sciImg1', 'sciImg2', 'sciImg3', 'sciImg4'])
    	#fftInfo_arg_df = fftInfo_arg.drop(['sciImg1', 'sciImg2', 'sciImg3', 'sciImg4'])
    	## ## for some reason, the dataframes make 3 identical rows; maybe because there are some (x,y) vectors
    	fftInfo_amp_df_this_frame = pd.DataFrame(fftInfo_amp)
    	fftInfo_arg_df_this_frame = pd.DataFrame(fftInfo_arg)
    	fftInfo_amp_df = fftInfo_amp_df.append(fftInfo_amp_df_this_frame, ignore_index=True)
    	fftInfo_arg_df = fftInfo_arg_df.append(fftInfo_arg_df_this_frame, ignore_index=True)

    	# print fyi
    	time_elapsed = time.time() - time_start
    	print("PSF "+file_name_full+" analyzed in time (secs):")
    	print(time_elapsed)

        # write to csv to check on a local machine (this info will also be pickled)
        fftInfo_amp_df.to_csv("pickled_info/csvs/fft_amp_"+str(int(counter_num))+".csv")
        fftInfo_arg_df.to_csv("pickled_info/csvs/fft_arg_"+str(int(counter_num))+".csv")
        ## ## note I havent used log_name anywhere yet

        ##################################################################################
        # write to csv in off-sky testing of retrieved FFT quantities using fake data
        csv_name = "trial1_190612_retrieved.csv"
        # find needed correction to FPC PL setpoint (degrees in K-band)
        #sci_to_K = np.divide(sci_wavel,2.2e-6) # factor to convert degrees in sci to degrees in K
        #corrxn_pl = -fftInfo_arg["med_highFreqPerfect_R"].values[0]*(180./np.pi)*sci_to_K
        #alpha_high_freq = [x_grad_perf_high_R, y_grad_perf_high_R] # gradient high freq lobe of PTF in x and y: [a,b]
        # state size of FFTed image
        Ny = np.shape(cookie_cut)[0]+2*padding_choice
        Nx = np.shape(cookie_cut)[1]+2*padding_choice
        # state gradient of PTF slope in x and y
        print("why median")
        print(fftInfo_arg["normVec_highFreqPerfect_R_x"])
        x_grad_perf_high_R = np.median(fftInfo_arg["normVec_highFreqPerfect_R_x"])
        y_grad_perf_high_R = np.median(fftInfo_arg["normVec_highFreqPerfect_R_y"])
        x_grad_perf_lowfreq = np.median(fftInfo_arg["normVec_lowFreqPerfect_x"])
        y_grad_perf_lowfreq = np.median(fftInfo_arg["normVec_lowFreqPerfect_y"])
        alpha_high_freq = [x_grad_perf_high_R, y_grad_perf_high_R] # gradient high freq lobe of PTF in x and y: [a,b]
        alpha_low_freq = [x_grad_perf_lowfreq, y_grad_perf_lowfreq] # same, in low freq lobe
        alpha_mean = np.mean([alpha_high_freq,alpha_low_freq],axis=0) # tip-tilt corrections should be based on gradients common to lobes (see Spalding+ SPIE 2018, Table 3)
        corrxn_tt = needed_tt_setpt_corrxn(alpha=alpha_mean,PS=plateScale_LMIR,Nx=Nx,Ny=Ny) # (x,y)
        with open(csv_name, "a") as datalist:
            opd_retrieve = fftInfo_arg["med_highFreqPerfect_R"]*((180./np.pi)/360.)*sci_wavel*1e6
            tip_retrieve = -corrxn_tt[1] # y (negative, because I want to know the retrieved value, not the correction)
            tilt_retrieve = -corrxn_tt[0] # x (ditto)
            # cols of macie time; fake file name; OPD; tip; tilt
            datalist.write("%f, %s, %f, %f, %f\n" % (time.time(), file_name_full, opd_retrieve, tip_retrieve, tilt_retrieve))

        # now, write out a
        # end code for writing retrieved values to csv
        ##################################################################################

        # pack info from the series of FFTs into dictionaries, and pickle them
        #d = {"fftInfo_amp": fftInfo_amp, "fftInfo_arg": fftInfo_arg}
        fft_pickle_write_name = "pickled_info/fft_info_"+str("{:0>2d}".format(counter_num))+".pkl"
        with open(fft_pickle_write_name, "w") as f:
            pickle.dump((fftInfo_amp, fftInfo_arg), f)

        print("Done analyzing one PSF.")
        counter_num += 1 # advance counter

    print("--------------------------")
    print("Done analyzing "+str(num_psfs_to_analyze)+" PSFs.")

    # return
    # 1. number of psfs over which we will take median
    # 2. the shape of the image we FFT so as to interpret the PTF slope
    return num_psfs_to_analyze, np.shape(amp)


def get_apply_pc_setpts(integ_time, num_psfs, fftimg_shape, sci_wavel, mode = "science", bkgd_mode = "quick_dirt", pickle_name = "setpts.pkl", apply = True):
    '''
    integ_time: vestigial, from when camera was being commanded by the Fizeau code
    num_psfs: number of PSFs to analyze
    fftimg_shape: shape of the array which actually gets FFTed
    sci_wavel: science wavelength
    mode: select testing or science mode
    apply=True: apply the corrections after they have been calculated; otherwise, just calculate them
    '''

    start_time = time.time()

    # restore the pickle file with the fit coefficients and scan data
    '''
    # for a single pickle file
    with open(fft_pickle_read_name) as f:
        fftInfo_amp, fftInfo_arg = pickle.load(f)
    #fftInfo_amp = d["fftInfo_amp"]
    #fftInfo_arg = d["fftInfo_arg"]
    '''

    # for multiple pickle files, over which we take medians of each key value
    d_all_amp = {}
    d_all_arg = {}
    for pickle_num in range(0,num_psfs): # loop over the number of PSFs which have been analyzed
        fft_pickle_read_name = "pickled_info/fft_info_"+str("{:0>2d}".format(pickle_num))+".pkl"

        # open individual pickle files (one per PSF) and put them into a larger dict
        # N.b. there are two dictionaries (for amp and arg) that come out of each pickle file
        with open(fft_pickle_read_name) as f:
            d_all_amp["pickle_" + str("{:0>2d}".format(pickle_num))], d_all_arg["pickle_" + str("{:0>2d}".format(pickle_num))] = pickle.load(f)

    # convert to one dataframe of dicts
    d_all_amp_df = pd.DataFrame(d_all_amp)
    d_all_arg_df = pd.DataFrame(d_all_arg)

    d_df_ampT = d_all_amp_df.T.reset_index() # transpose
    d_df_argT = d_all_arg_df.T.reset_index() # transpose
    del d_df_ampT["index"], d_df_ampT["dummy"] # get rid of dummy so I can take a median
    del d_df_argT["index"], d_df_argT["dummy"] # get rid of dummy so I can take a median

    print('----------------------------')
    # find median values of quantities across the dicts
    fftInfo_amp = pd.DataFrame(pd.DataFrame.median(d_df_ampT, axis=0)).T
    fftInfo_arg = pd.DataFrame(pd.DataFrame.median(d_df_argT, axis=0)).T
    #fftInfo_amp = pd.DataFrame(pd.DataFrame.median(d_all_amp_df, axis=1, skipna=True)).T
    #fftInfo_arg = pd.DataFrame(pd.DataFrame.median(d_all_arg_df, axis=1, skipna=True)).T
    print(fftInfo_amp)
    print(fftInfo_arg)
    #fftInfo_amp = d["fftInfo_amp"]
    #fftInfo_arg = d["fftInfo_arg"]

    #####################################
    # LOW PRIORITY: Overlap of the Airy PSFs (MAYBE THIS SHOULD BE RESOLVED WITH OVERLAP_PSFS?)
    # TO CORRECT: TT THE HPC TO TRANSLATE DX AIRY PSF
    ## ## (MAYBE TRY TTING THE FPC LATER?)
    ## MAYBE TRY MOVING STUFF MANUALLY ON PIEZOS PAGE FIRST
    print("----------------------------------")
    print("Checking overlap of the Airy PSFs via std of FFT phase low freq node:")
    print(np.median(fftInfo_arg["std_lowFreqPerfect"]))
    #if (mode != "total_passive"):
    if True: # total_passive is not a mode now
        hpc_tip_position_now = pi.getINDI("dac_stage.hpc_status.tip") # HPC tip pos
        hpc_tilt_position_now = pi.getINDI("dac_stage.hpc_status.tilt") # HPC tilt pos
        fpc_tip_setpoint_now = pi.getINDI("dac_stage.fpc_status.tip") # FPC tip setpoint
        fpc_tilt_setpoint_now = pi.getINDI("dac_stage.fpc_status.tilt") # FPC tilt setpoint
        print("----------------------------------")
        print("hpc_tip_position_now:")
        print(hpc_tip_position_now)
        print("hpc_tilt_position_now:")
        print(hpc_tilt_position_now)
        print("fpc_tip_setpoint_now:")
        print(fpc_tip_setpoint_now)
        print("fpc_tilt_setpoint_now:")
        print(fpc_tilt_setpoint_now)
        print("Check the Airy overlap on LMIR, and manually move the HPC. How large of an HPC tip-tilt was necessary to overlap them? Scale stdev with this movement.")

    #####################################
    # MEDIUM PRIORITY: High-freq fringe visibility (as measued with median of high-freq FFT ampl lobe) (ASSUME OPD=0 BEFORE TAKING OUT OTHER ABERRATIONS)
    print("----------------------------------")
    print("Checking high-freq fringe visibility via median of ampl of high freq lobe in MTF:")
    print(np.median(fftInfo_amp["med_highFreqPerfect_R"]))
    print("----------------------------------")
    #if (mode != "total_passive"):
    if True: # total_passive is not a mode now
        print("Current SPC_Trans position:")
        spc_trans_position_now = pi.getINDI("Ubcs.SPC_Trans_status.PosNum") # translation stage (absolute position, 0.02 um)
        print(spc_trans_position_now)
        print("Current FPC piston position:")
        print("Current FPC PL setpoint:")
        fpc_pl_setpoint = pi.getINDI("PLC.UBCSettings.PLSetpoint")
        print("Manually change PL setpoint. How good/bad is the high-freq fringe visibility on LMIR? Write down the scale.")

    #####################################
    # HIGH PRIORITY: High-freq phase gradient
    print("----------------------------------")
    print("Checking ampl of high-freq node of PTF:")
    print(fftInfo_arg["med_highFreqPerfect_R"].values[0])
    print("Checking phase gradient in x of high freq in PTF:")
    x_grad_perf_high_R = np.median(fftInfo_arg["normVec_highFreqPerfect_R_x"])
    print(x_grad_perf_high_R)
    print("Checking phase gradient in y of high freq in PTF:")
    y_grad_perf_high_R = np.median(fftInfo_arg["normVec_highFreqPerfect_R_y"])
    print(y_grad_perf_high_R)
    print("Checking phase gradient in x of low freq in PTF:")
    x_grad_perf_lowfreq = np.median(fftInfo_arg["normVec_lowFreqPerfect_x"])
    print(x_grad_perf_lowfreq)
    print("Checking phase gradient in y of low freq in PTF:")
    y_grad_perf_lowfreq = np.median(fftInfo_arg["normVec_lowFreqPerfect_y"])
    print(y_grad_perf_lowfreq)
    print("Checking phase gradient in x,y (amplitude) of high freq in PTF:")
    print(np.sqrt(np.add(math.pow(x_grad_perf_high_R,2),math.pow(y_grad_perf_high_R,2))))
    print("Checking angle of the gradient (deg CCW from x=0) in PTF:")
    angle_gradient = math.atan2(y_grad_perf_high_R,x_grad_perf_high_R)*180./np.pi
    print(angle_gradient)
    print("--------------------------")

    # find needed correction to FPC TT setpoints (mas)
    alpha_high_freq = [x_grad_perf_high_R, y_grad_perf_high_R] # gradient high freq lobe of PTF in x and y: [a,b]
    print('alpha_high_freq:')
    print(alpha_high_freq)
    alpha_low_freq = [x_grad_perf_lowfreq, y_grad_perf_lowfreq] # same, in low freq lobe
    print('alpha_low_freq:')
    print(alpha_low_freq)
    alpha_mean = np.mean([alpha_high_freq,alpha_low_freq],axis=0) # corrections should be based on gradients common to lobes (see Spalding+ SPIE 2018, Table 3)
    print('alpha mean:')
    print(alpha_mean) # alpha_mean = [alpha_x from both low- and high-freq nodes, alpha_y from both low- and high-freq nodes]
    Nx = fftimg_shape[1]
    Ny = fftimg_shape[0]
    import pdb; pdb.set_trace()
    corrxn_tt = needed_tt_setpt_corrxn(alpha=alpha_mean,PS=plateScale_LMIR,Nx=Nx,Ny=Ny) # (x,y)
    corrxn_tilt_x = int(corrxn_tt[0])
    corrxn_tip_y = int(corrxn_tt[1])
    print("Needed TT setpt correction (x, y) in mas")
    print(corrxn_tt)

    # find needed correction to FPC PL setpoint (degrees in K-band)
    # (even if there is a slope in the PTF, the median would be zero if PL error is zero)
    sci_to_K = np.divide(sci_wavel,2.2e-6) # factor to convert degrees in sci to degrees in K
    corrxn_pl = -fftInfo_arg["med_highFreqPerfect_R"].values[0]*(180./np.pi)*sci_to_K

    print("Current FPC tip setpoint:")
    fpc_tip_setpoint = pi.getINDI("PLC.UBCSettings.TipSetpoint")
    print(fpc_tip_setpoint)
    print("Current FPC tilt setpoint:")
    fpc_tilt_setpoint = pi.getINDI("PLC.UBCSettings.TiltSetpoint")
    print(fpc_tilt_setpoint)
    print("Current FPC pl setpoint:")
    fpc_pl_setpoint = pi.getINDI("PLC.UBCSettings.PLSetpoint")
    print(fpc_pl_setpoint)
    print("-----")
    print("FPC tip (y) setpoint correction (mas):")
    print(str(corrxn_tip_y))
    print("FPC tilt (x) setpoint correction (mas):")
    print(str(corrxn_tilt_x))
    print("FPC pl (piston) setpoint correction (deg):")
    print(str(corrxn_pl))
    print("-----")

    # get other values that need to be included in the setINDI command
    PLC_UBCSettings_Beam1_x = pi.getINDI("PLC.UBCSettings.Beam1_x")
    PLC_UBCSettings_Beam1_y = pi.getINDI("PLC.UBCSettings.Beam1_y")
    PLC_UBCSettings_Beam2_x = pi.getINDI("PLC.UBCSettings.Beam2_x")
    PLC_UBCSettings_Beam2_y = pi.getINDI("PLC.UBCSettings.Beam2_y")
    PLC_UBCSettings_Beam_r = pi.getINDI("PLC.UBCSettings.Beam_r")
    PLC_UBCSettings_MinFTSNR = pi.getINDI("PLC.UBCSettings.MinFTSNR")
    PLC_UBCSettings_PWVGain = pi.getINDI("PLC.UBCSettings.PWVGain")
    PLC_UBCSettings_PLPGain = pi.getINDI("PLC.UBCSettings.PLPGain")
    PLC_UBCSettings_PLDGain = pi.getINDI("PLC.UBCSettings.PLDGain")
    PLC_UBCSettings_PLIGain = pi.getINDI("PLC.UBCSettings.PLIGain")
    PLC_UBCSettings_CGSetpoint = pi.getINDI("PLC.UBCSettings.CGSetpoint")
    PLC_UBCSettings_CGScale = pi.getINDI("PLC.UBCSettings.CGScale")
    PLC_UBCSettings_TTPGain = pi.getINDI("PLC.UBCSettings.TTPGain")
    PLC_UBCSettings_TTDGain = pi.getINDI("PLC.UBCSettings.TTDGain")
    PLC_UBCSettings_TTIGain = pi.getINDI("PLC.UBCSettings.TTIGain")
    ''' UNCOMMENT HERE
    if (apply == True):
        # all these have to be applied at once; note that tip is y, tilt is x
        # print("Manually change FPC tip-tilt setpoints. What was the scale?")
        pi.setINDI("PLC.UBCSettings.Beam1_x="+str(PLC_UBCSettings_Beam1_x)+";Beam1_y="+str(PLC_UBCSettings_Beam1_y)+";Beam2_x="+str(PLC_UBCSettings_Beam2_x)+";Beam2_y="+str(PLC_UBCSettings_Beam2_y)+";Beam_r="+str(PLC_UBCSettings_Beam_r)+";MinFTSNR="+str(PLC_UBCSettings_MinFTSNR)+";PLSetpoint="+str(np.add(fpc_pl_setpoint,corrxn_pl)))+";PWVGain="+str(PLC_UBCSettings_PWVGain)+";PLIGain="+str(PLC_UBCSettings_PLIGain)+";CGSetpoint="+str(PLC_UBCSettings_CGSetpoint)+";CGScale="+str(PLC_UBCSettings_CGScale)+";TipSetpoint="+str(np.add(fpc_tip_setpoint,corrxn_tip_y))+";TiltSetpoint="+str(np.add(fpc_tilt_setpoint,corrxn_tilt_x))+";TTIGain="+str(PLC_UBCSettings_TTIGain)+";PLPGain="+str(PLC_UBCSettings_PLPGain)+";PLDGain="+str(PLC_UBCSettings_PLDGain)+";TTPGain="+str(PLC_UBCSettings_TTPGain)+";TTDGain="+str(PLC_UBCSettings_TTDGain))
    '''
    #######################################################################
    # lines to run on the command line to test application of corrections
    '''
    # LOW PRIORITY: To correct Airy PSF overlap, or TT in Fizeau PSF
    new_tip_setpoint = 0
    new_tilt_setpoint = 0
    if ((mode == "science") or (mode == "fake_fits") or (mode == "az_source")):
        pi.setINDI("PLC.UBCSettings.Beam1_x="+str(PLC_UBCSettings_Beam1_x)+";Beam1_y="+str(PLC_UBCSettings_Beam1_y)+";Beam2_x="+str(PLC_UBCSettings_Beam2_x)+";Beam2_y="+str(PLC_UBCSettings_Beam2_y)+";Beam_r="+str(PLC_UBCSettings_Beam_r)+";MinFTSNR="+str(PLC_UBCSettings_MinFTSNR)+";PLSetpoint="+str(fpc_pl_setpoint)+";PWVGain="+str(PLC_UBCSettings_PWVGain)+";PLIGain="+str(PLC_UBCSettings_PLIGain)+";CGSetpoint="+str(PLC_UBCSettings_CGSetpoint)+";CGScale="+str(PLC_UBCSettings_CGScale)+";TipSetpoint="+str(int(new_tip_setpoint))+";TiltSetpoint="+str(int(new_tilt_setpoint))+";TTIGain="+str(PLC_UBCSettings_TTIGain)+";PLPGain="+str(PLC_UBCSettings_PLPGain)+";PLDGain="+str(PLC_UBCSettings_PLDGain)+";TTPGain="+str(PLC_UBCSettings_TTPGain)+";TTDGain="+str(PLC_UBCSettings_TTDGain))

    # MEDIUM PRIORITY: To correct high-freq fringe visibility
    new_pl_setpoint = 0
    spc_trans_stepSize = 5. # (um, total OPD)
    if ((mode == "science") or (mode == "fake_fits") or (mode == "az_source")):
        pi.setINDI("PLC.UBCSettings.PLSetpoint="+str(int(new_pl_setpoint)))
        pi.setINDI("PLC.UBCSettings.Beam1_x="+str(PLC_UBCSettings_Beam1_x)+";Beam1_y="+str(PLC_UBCSettings_Beam1_y)+";Beam2_x="+str(PLC_UBCSettings_Beam2_x)+";Beam2_y="+str(PLC_UBCSettings_Beam2_y)+";Beam_r="+str(PLC_UBCSettings_Beam_r)+";MinFTSNR="+str(PLC_UBCSettings_MinFTSNR)+";PLSetpoint="+str(int(new_pl_setpoint))+";PWVGain="+str(PLC_UBCSettings_PWVGain)+";PLIGain="+str(PLC_UBCSettings_PLIGain)+";CGSetpoint="+str(PLC_UBCSettings_CGSetpoint)+";CGScale="+str(PLC_UBCSettings_CGScale)+";TipSetpoint="+str(fpc_tip_setpoint)+";TiltSetpoint="+str(fpc_tilt_setpoint)+";TTIGain="+str(PLC_UBCSettings_TTIGain)+";PLPGain="+str(PLC_UBCSettings_PLPGain)+";PLDGain="+str(PLC_UBCSettings_PLDGain)+";TTPGain="+str(PLC_UBCSettings_TTPGain)+";TTDGain="+str(PLC_UBCSettings_TTDGain))
        pi.setINDI("Ubcs.SPC_Trans.command=>"+'{0:.1f}'.format(spc_trans_command))
        pi.setINDI("Ubcs.SPC_Trans.command=>"+'{0:.1f}'.format(10*0.5*stepSize)) # factor of 10 bcz command is in 0.1 um

    # HIGH PRIORITY: To correct high-freq fringe gradients (note similarity to correction for Airy PSF overlap)
    new_tip_setpoint = 0
    new_tilt_setpoint = 0
    if ((mode == "science") or (mode == "fake_fits") or (mode == "az_source")):
        pi.setINDI("PLC.UBCSettings.TipSetpoint="+str(int(new_tip_setpoint)))
        pi.setINDI("PLC.UBCSettings.TiltSetpoint="+str(int(new_tilt_setpoint)))

    # FYI: EDIT TO MAKE MOMENTARY CHANGES IN FPC TT
    if ((mode == "science") or (mode == "fake_fits") or (mode == "az_source")):
        pi.setINDI("dac_stage.fpc.tip="+'{0:.1f}'.format(vector_move_asec[0])+";tilt="+'{0:.1f}'.format(vector_move_asec[1])+";piston=0;mode=1")
    '''

    '''UNCOMMENT HERE
    # pickle the calculated corrections, so that they can be checked in the next function
    with open(pickle_name, "w") as f:
        pickle.dump((corrxn_tt,corrxn_pl), f)

    # print needed corrections
    print("--------------------------------------")
    print("--------------------------------------")
    print("NEEDED TIP (y) CORRECTION TO FPC SETPOINT (mas):")
    print(str(int(corrxn_tt[1])))
    print("---")
    print("NEEDED TILT (x) CORRECTION TO FPC SETPOINT (mas):")
    print(str(int(corrxn_tt[0])))
    print("---")
    print("NEEDED PL CORRECTION TO FPC SETPOINT (degrees):")
    # positive number means jailbars should be moved to the right on LMIR by this amount
    print(str(corrxn_pl))

    stop_time = time.time()
    print("-")
    print("Elapsed time for sending setpoint or PL corrections:")
    print(stop_time - start_time)
    print("-----------------------------------")

    # turn off fizeau flag to avoid problems with other observations
    print("De-activating ROI aquisition flag")
    pi_fiz.setINDI("fizeau.enable_run.value=Off")
    '''
    return


def compare_setpts(pickle_pre, pickle_post, mode = "science", bkgd_mode = "quick_dirt"):
    '''
    Check that the PSF has been improved with the new setpoints. (Due to sign degeneracies,
    they might have made things worse after the first try.)

    INPUTS:
    pickle_pre: filename of pickle file containing TT,PL setpoint corrections after the 1st correction
    pickle_post: same as above, except its for the needed corrections AFTER the 1st was made

    OUTPUTS:
    '''

    # restore the pickle file containing the 1st calculated corrections
    with open(pickle_pre) as f:
        corrxn_tt_pre, corrxn_pl_pre = pickle.load(f)

    # restore pickle file with needed corrections after the 1st correction
    with open(pickle_post) as g:
        corrxn_tt_post, corrxn_pl_post = pickle.load(g)

    # print info
    print("-----------------")
    print("Needed tip (y) correction, pre and post the application of the first correction:")
    print("Pre: " + str(corrxn_tt_pre[0]))
    print("Post: " + str(corrxn_tt_post[0]))
    print("Needed tilt (x) correction, pre and post the application of the first correction:")
    print("Pre: " + str(corrxn_tt_pre[1]))
    print("Post: " + str(corrxn_tt_post[1]))
    print("Needed pathlength correction, pre and post the application of the first correction:")
    print("Pre: " + str(corrxn_pl_pre))
    print("Post: " + str(corrxn_pl_post))

    # get other values that need to be included in the setINDI command
    PLC_UBCSettings_Beam1_x = pi.getINDI("PLC.UBCSettings.Beam1_x")
    PLC_UBCSettings_Beam1_y = pi.getINDI("PLC.UBCSettings.Beam1_y")
    PLC_UBCSettings_Beam2_x = pi.getINDI("PLC.UBCSettings.Beam2_x")
    PLC_UBCSettings_Beam2_y = pi.getINDI("PLC.UBCSettings.Beam2_y")
    PLC_UBCSettings_Beam_r = pi.getINDI("PLC.UBCSettings.Beam_r")
    PLC_UBCSettings_MinFTSNR = pi.getINDI("PLC.UBCSettings.MinFTSNR")
    PLC_UBCSettings_PWVGain = pi.getINDI("PLC.UBCSettings.PWVGain")
    PLC_UBCSettings_PLPGain = pi.getINDI("PLC.UBCSettings.PLPGain")
    PLC_UBCSettings_PLDGain = pi.getINDI("PLC.UBCSettings.PLDGain")
    PLC_UBCSettings_PLIGain = pi.getINDI("PLC.UBCSettings.PLIGain")
    PLC_UBCSettings_CGSetpoint = pi.getINDI("PLC.UBCSettings.CGSetpoint")
    PLC_UBCSettings_CGScale = pi.getINDI("PLC.UBCSettings.CGScale")
    PLC_UBCSettings_TTPGain = pi.getINDI("PLC.UBCSettings.TTPGain")
    PLC_UBCSettings_TTDGain = pi.getINDI("PLC.UBCSettings.TTDGain")
    PLC_UBCSettings_TTIGain = pi.getINDI("PLC.UBCSettings.TTIGain")

    # all these have to be applied at once; note that tip is y, tilt is x
    '''UNCOMMENT HERE
    pi.setINDI("PLC.UBCSettings.Beam1_x="+str(PLC_UBCSettings_Beam1_x)+";Beam1_y="+str(PLC_UBCSettings_Beam1_y)+";Beam2_x="+str(PLC_UBCSettings_Beam2_x)+";Beam2_y="+str(PLC_UBCSettings_Beam2_y)+";Beam_r="+str(PLC_UBCSettings_Beam_r)+";MinFTSNR="+str(PLC_UBCSettings_MinFTSNR)+";PLSetpoint="+str(np.add(fpc_pl_setpoint,corrxn_pl)))+";PWVGain="+str(PLC_UBCSettings_PWVGain)+";PLIGain="+str(PLC_UBCSettings_PLIGain)+";CGSetpoint="+str(PLC_UBCSettings_CGSetpoint)+";CGScale="+str(PLC_UBCSettings_CGScale)+";TipSetpoint="+str(np.add(fpc_tip_setpoint,corrxn_tip_y))+";TiltSetpoint="+str(np.add(fpc_tilt_setpoint,corrxn_tilt_x))+";TTIGain="+str(PLC_UBCSettings_TTIGain)+";PLPGain="+str(PLC_UBCSettings_PLPGain)+";PLDGain="+str(PLC_UBCSettings_PLDGain)+";TTPGain="+str(PLC_UBCSettings_TTPGain)+";TTDGain="+str(PLC_UBCSettings_TTDGain))

    # if not, make a 2x reverse correction
    # check if correction has gone further away from zero
    if ((  np.abs(corrxn_tt_post[0]) > np.abs(corrxn_tt_pre[0])  ) and (  np.sign(corrxn_tt_post[0]) == np.sign(corrxn_tt_pre[0])  )):
        print("Re-correcting the tip (y) setpoint correction")
        fpc_tip_setpoint = pi.getINDI("PLC.UBCSettings.TipSetpoint")
        fpc_tilt_setpoint = pi.getINDI("PLC.UBCSettings.TiltSetpoint")
        fpc_pl_setpoint = pi.getINDI("PLC.UBCSettings.PLSetpoint")
        new_tip_offset = -2*corrxn_tt_pre[0]
        pi.setINDI("PLC.UBCSettings.Beam1_x="+str(PLC_UBCSettings_Beam1_x)+";Beam1_y="+str(PLC_UBCSettings_Beam1_y)+";Beam2_x="+str(PLC_UBCSettings_Beam2_x)+";Beam2_y="+str(PLC_UBCSettings_Beam2_y)+";Beam_r="+str(PLC_UBCSettings_Beam_r)+";MinFTSNR="+str(PLC_UBCSettings_MinFTSNR)+";PLSetpoint="+str(fpc_pl_setpoint)+";PWVGain="+str(PLC_UBCSettings_PWVGain)+";PLIGain="+str(PLC_UBCSettings_PLIGain)+";CGSetpoint="+str(PLC_UBCSettings_CGSetpoint)+";CGScale="+str(PLC_UBCSettings_CGScale)+";TipSetpoint="+str(np.add(fpc_tip_setpoint,new_tip_offset))+";TiltSetpoint="+str(fpc_tilt_setpoint)+";TTIGain="+str(PLC_UBCSettings_TTIGain)+";PLPGain="+str(PLC_UBCSettings_PLPGain)+";PLDGain="+str(PLC_UBCSettings_PLDGain)+";TTPGain="+str(PLC_UBCSettings_TTPGain)+";TTDGain="+str(PLC_UBCSettings_TTDGain))
    if ((  np.abs(corrxn_tt_post[1]) > np.abs(corrxn_tt_pre[1])  ) and (  np.sign(corrxn_tt_post[1]) == np.sign(corrxn_tt_pre[1])  )):
        print("Re-correcting the tilt (x) setpoint correction")
        fpc_tip_setpoint = pi.getINDI("PLC.UBCSettings.TipSetpoint")
        fpc_tilt_setpoint = pi.getINDI("PLC.UBCSettings.TiltSetpoint")
        fpc_pl_setpoint = pi.getINDI("PLC.UBCSettings.PLSetpoint")
        new_tilt_offset = -2*corrxn_tt_pre[1]
        pi.setINDI("PLC.UBCSettings.Beam1_x="+str(PLC_UBCSettings_Beam1_x)+";Beam1_y="+str(PLC_UBCSettings_Beam1_y)+";Beam2_x="+str(PLC_UBCSettings_Beam2_x)+";Beam2_y="+str(PLC_UBCSettings_Beam2_y)+";Beam_r="+str(PLC_UBCSettings_Beam_r)+";MinFTSNR="+str(PLC_UBCSettings_MinFTSNR)+";PLSetpoint="+str(fpc_pl_setpoint)+";PWVGain="+str(PLC_UBCSettings_PWVGain)+";PLIGain="+str(PLC_UBCSettings_PLIGain)+";CGSetpoint="+str(PLC_UBCSettings_CGSetpoint)+";CGScale="+str(PLC_UBCSettings_CGScale)+";TipSetpoint="+str(fpc_tip_setpoint)+";TiltSetpoint="+str(np.add(fpc_tilt_setpoint,new_tilt_offset))+";TTIGain="+str(PLC_UBCSettings_TTIGain)+";PLPGain="+str(PLC_UBCSettings_PLPGain)+";PLDGain="+str(PLC_UBCSettings_PLDGain)+";TTPGain="+str(PLC_UBCSettings_TTPGain)+";TTDGain="+str(PLC_UBCSettings_TTDGain))
    if ((  np.abs(corrxn_pl_post) > np.abs(corrxn_pl_pre)  ) and (  np.sign(corrxn_pl_post) == np.sign(corrxn_pl_pre)  )):
        print("Re-correcting the PL setpoint correction")
        fpc_tip_setpoint = pi.getINDI("PLC.UBCSettings.TipSetpoint")
        fpc_tilt_setpoint = pi.getINDI("PLC.UBCSettings.TiltSetpoint")
        fpc_pl_setpoint = pi.getINDI("PLC.UBCSettings.PLSetpoint")
        new_pl_offset = -2*corrxn_pl_pre
        pi.setINDI("PLC.UBCSettings.Beam1_x="+str(PLC_UBCSettings_Beam1_x)+";Beam1_y="+str(PLC_UBCSettings_Beam1_y)+";Beam2_x="+str(PLC_UBCSettings_Beam2_x)+";Beam2_y="+str(PLC_UBCSettings_Beam2_y)+";Beam_r="+str(PLC_UBCSettings_Beam_r)+";MinFTSNR="+str(PLC_UBCSettings_MinFTSNR)+";PLSetpoint="+str(np.add(fpc_pl_setpoint,new_pl_offset))+";PWVGain="+str(PLC_UBCSettings_PWVGain)+";PLIGain="+str(PLC_UBCSettings_PLIGain)+";CGSetpoint="+str(PLC_UBCSettings_CGSetpoint)+";CGScale="+str(PLC_UBCSettings_CGScale)+";TipSetpoint="+str(fpc_tip_setpoint)+";TiltSetpoint="+str(fpc_tilt_setpoint)+";TTIGain="+str(PLC_UBCSettings_TTIGain)+";PLPGain="+str(PLC_UBCSettings_PLPGain)+";PLDGain="+str(PLC_UBCSettings_PLDGain)+";TTPGain="+str(PLC_UBCSettings_TTPGain)+";TTDGain="+str(PLC_UBCSettings_TTDGain))
    '''
