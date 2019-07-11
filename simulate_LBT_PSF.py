#from pathos.multiprocessing import ProcessingPool as Pool

# coding: utf-8
import pickle
#import jFits
import numpy as np
import matplotlib.pyplot as mpl
from astropy.io import fits
from os import listdir
import glob
import multiprocessing as mp
from multiprocessing import Process, Queue, Pool
import itertools
import scipy
#from scipy import ndimage
import scipy.ndimage
import photutils
import time
import astropy
from astropy.io import fits
import numpy.ma as ma
import matplotlib.pyplot as plt
from walker_def import walker


def lbt_aperture2(deltaspacing,N_pix,D):
    # create LBT aperture
    #deltaspacing = 0.100 # 62.5 #mm; spacing of pixels in the aperture
    #N_pix = 2048 # samples across the aperture image
    #D = 8.25 # M1 diameter, in m
    padfactor = N_pix*deltaspacing/D
    
    centerx=(D*padfactor/deltaspacing)/2.0 #pixels
    centery=(D*padfactor/deltaspacing)/2.0 #pixels
    lenx = D*padfactor/deltaspacing #pixels
    leny = D*padfactor/deltaspacing #pixels
    hwstamp = D/2.0/deltaspacing #pixels

    xpos = np.arange(-D/2.0*padfactor,D/2.0*padfactor,deltaspacing)
    ypos = np.arange(-D/2.0*padfactor,D/2.0*padfactor,deltaspacing)
    xx,yy = np.meshgrid(xpos,ypos)

    # SX, DX centers
    xsx = xx - 7.2085
    xdx = xx + 7.2085

    # define apertures within radii
    rr=np.sqrt(xx**2+yy**2)
    rrsx= np.sqrt(xsx**2+yy**2)
    rrdx= np.sqrt(xdx**2+yy**2)
    sxAperture=rrsx<D/2.0
    dxAperture=rrdx<D/2.0
   
    LBTAperture = sxAperture + dxAperture
    
    return centerx, centery, sxAperture, dxAperture, LBTAperture, deltaspacing

###################################


def make_psf(sxAperture, dxAperture, phase, wavel, centerx, centery, tip_asec, tilt_asec):
    # make (complex) psf from aperture and phasor

    # aperture plate scale has to be re-computed for every wavelength
    aperturePlateScale = 9416.20*wavel
    
    # make phase screen for injecting relative tip-tilt
    phaseScreenIndicesUnnormzed = np.indices(np.shape(sxAperture)) # use indices to make screen
    # if the tip (or tilt) value in the aperture center is zero, what value
    tip_unnormzed = aperturePlateScale*phaseScreenIndicesUnnormzed[0]*tip_asec/206264.806 # small angle approx for np.tan(tip_asec/206264.806)
    tilt_unnormzed = aperturePlateScale*phaseScreenIndicesUnnormzed[1]*tilt_asec/206264.806 
    
    tipDXavg = np.mean(tip_unnormzed[np.where(dxAperture)]) #  mean value of the tip/tilt in the DX aperture
    tiltDXavg = np.mean(tilt_unnormzed[np.where(dxAperture)])
    
    tip = np.divide(np.subtract(tip_unnormzed,tipDXavg),wavel) # make the OPD from tip equal zero at the center in the DX aperture, and don't forget to divide by lambda
    tilt = np.divide(np.subtract(tilt_unnormzed,tiltDXavg),wavel) # make piston-neutral and div by lambda

    phaseScreen = np.add(tip,tilt) # adds x- and y- components of screen, i.e., tip and tilt
    
    sxAperturephase = np.multiply(1.0, sxAperture)
    dxAperturephase = np.multiply(dxAperture, np.exp(1.j * 2 * np.pi * np.add(phase,phaseScreen))) # phase: for OPD; phaseScreen: for tip-tilt
    
    LBTAperturephase = sxAperturephase + dxAperturephase
    
    FTAp = np.fft.fft2(LBTAperturephase)
    FTAp = np.fft.fftshift(FTAp)

    return FTAp


###################################


def psf_sim(parameters):
    # opd in meters
    # diff_tip and diff_tilt in asec
    # transl in pixels

    
    start_time = time.time()

    #################################
    ####### BEGIN USER INPUTS #######
    
    monochromatic = True
    lmir = False # lmircam or nomic?
    
    if monochromatic:
        wavel = 12.0e-6 # m (monochromatic)
    else:
        wavel = np.linspace(4.02,4.08,num=10)*1e-6 # (polychromatic)
        
    ####### END USER INPUTS #######
    #################################
        
    if lmir:
        omega = 9412.639 # LMIR plate scale is 0.0107 asec/pix (see research journal from 2018 Apr 19 or 2019 July 11)
        plateScale = 0.0107 # in asec/pix
    else:
        omega = 5595.291 # NOMIC plate scale is 0.018 asec/pix
        plateScale = 0.018 # in asec/pix
    
    
    N_pix = 2048 # number of pixels across the input and output arrays
    D = 8.25 # effective M1 diameter, in m
    avg_wavel = np.mean(wavel)
    extraString = '' # some string to distinguish these pickle files

    opd = parameters[0]
    diff_tip = parameters[1]
    diff_tilt = parameters[2]
    transl = parameters[3]
    index_num = parameters[4]
    
    # dictionary of psf stats
    psf_stats = {'OPD': [], 'tip': [], 'tilt': [], 'translation': [],
             'centerlobe_ampl': [], 'sidelobe_ampl': [],
             'centerlobe_phase_deg': [], 'sidelobe_phase_deg': [],
             'PSF_image': [], 'FTamp_image': [], 'FTphase_image': [],
             'freq_x_axis_shifted': [], 'freq_y_axis_shifted': []} 
    
    phase = opd/wavel # this is a float (monochrom) or array (polychrom)
    
    # INJECT RELATIVE TIP-TILT
    tip_asec = diff_tip
    tilt_asec = diff_tilt
        
    # create mono- or polychromatic PSF
    if monochromatic: # wavel just has 1 value
        deltaSpacing = omega*wavel # sampling rate in pupil plane (m/pix)
        centerx, centery, sxAperture, dxAperture, LBTAperture, aperturePlateScale = lbt_aperture2(deltaSpacing, N_pix, D) # need to re-generate LBT aperture for each wavelength for DFT
        FTAp = make_psf(sxAperture, dxAperture, phase, wavel, centerx, centery, tip_asec, tilt_asec)
        I = np.real(FTAp * np.conj(FTAp))
    else: # polychromatic PSF
        for t in range(0,len(wavel)): # integrate over wavelength
            print('Polychromatic PSF, calculating phase = '+str(phase[t]))
            deltaSpacing = omega*wavel # sampling rate in pupil plane (m/pix)
            centerx, centery, sxAperture, dxAperture, LBTAperture, aperturePlateScale = lbt_aperture2(deltaSpacing[t], N_pix, D) # need to re-generate LBT aperture for each wavelength for DFT
            FTAp = make_psf(sxAperture, dxAperture, phase[t], wavel[t], centerx, centery, tip_asec, tilt_asec)
            if (t==0):
                Iall = np.zeros(np.shape(FTAp)) # initial array of zeros
            I = np.real(FTAp * np.conj(FTAp))
            Iall = np.add(Iall,I) # pile illumination from different wavelengths on top of each other at the center (the FT images will be of different sizes)
        I = np.copy(Iall)

    # this needed to get rid of the checkerboard pattern
    cutI = I # (vestigial)
    padI = cutI # np.pad(cutI,len(cutI),'edge') for addl padding
    padI_no_fft_shift_0 = padI # needed for plotting later
    padI = np.fft.fftshift(padI)

    
    # INJECT PSF TRANSLATION
    #if any_y_translation: # fixed y-translation
    #    transl_y = 1.0
    #else:
    #    transl_y = 0.0
    transl_y = 0.0 # zero for now
    padI_no_translation = np.copy(padI)
    padI = scipy.ndimage.interpolation.shift(padI_no_translation,[transl_y,transl])
    padI_transl_but_no_fft_shift = scipy.ndimage.interpolation.shift(padI_no_fft_shift_0,[transl_y,transl])

    # print stuff
    print('----------------')
    print('Making PSF for ')
    print('OPD: '+str(opd*1e6)+' um')
    print('Tip: '+str(tip_asec)+' asec')
    print('Tilt: '+str(tilt_asec)+' asec')
    print('Translation in x: '+str(transl)+' pix')
    print('Translation in y: '+str(transl_y)+' pix')
                
    # now extract results from the PSF...
                
    PhaseExtract = np.fft.fft2(padI)
    PhaseExtract = np.fft.fftshift(PhaseExtract)
    
    AmpPE = np.absolute(PhaseExtract)
    ArgPE = np.angle(PhaseExtract)
    ArgPE_deg = ArgPE*180./np.pi
    CenterPix = len(PhaseExtract)/2.0
    IntPix = 150.0;

    # find amplitude of
    # central lobe of FT amplitude
    center_ampl = np.max(ma.masked_array(AmpPE, mask=mask_center))
    # side lobe of FT amplitude
    side_ampl = np.max(ma.masked_array(AmpPE, mask=mask_right))
    
    # find phase value at center and side lobe
    
    # ... central lobe of FT amplitude
    CenterPhase_deg = np.ma.median(ma.masked_array(ArgPE_deg, mask=mask_center))
    sidelobe_phase_deg = np.ma.median(ma.masked_array(ArgPE_deg, mask=mask_right))
    
    # get the axes right
    freq_x_axis = np.fft.fftfreq(np.shape(I)[1]) # frequency axes (in pix^-1)
    freq_y_axis = np.fft.fftfreq(np.shape(I)[0])
    freq_x_axis_shifted = np.fft.fftshift(freq_x_axis) 
    freq_y_axis_shifted = np.fft.fftshift(freq_y_axis)
    
    # store stats and images
    psf_stats['OPD'] = opd
    psf_stats['translation'] = transl
    psf_stats['tip'] = tip_asec
    psf_stats['tilt'] = tilt_asec
    psf_stats['sidelobe_phase_deg'] = sidelobe_phase_deg     
    psf_stats['centerlobe_phase_deg'] = CenterPhase_deg
    psf_stats['PSF_image'] = padI_transl_but_no_fft_shift
    psf_stats['FTamp_image'] = AmpPE
    psf_stats['FTphase_image'] = ArgPE_deg
    psf_stats['centerlobe_ampl'] = center_ampl
    psf_stats['sidelobe_ampl'] = side_ampl
    psf_stats['freq_x_axis_shifted'] = freq_x_axis_shifted
    psf_stats['freq_y_axis_shifted'] = freq_y_axis_shifted                 

    # strings for filenames
    if monochromatic:
        chromaticString = 'mono'
    else:
        chromaticString = 'poly'
    wavelString = str(int(avg_wavel*1e9))
    opdString = str(int(opd*1e8)).zfill(5)
    tipString = str(int(1000*diff_tip)).zfill(4)
    tiltString = str(int(1000*diff_tilt)).zfill(4)
    translString = str(int(100*transl)).zfill(3)
    PSstring = str(int(1000*plateScale))

    #mpl.imshow(ArgPE_deg) # just FYI
    #mpl.savefig('test.png')

    #hdu = fits.PrimaryHDU(ArgPE_deg) # just FYI
    #hdu.writeto('test.fits', overwrite=True)

    # save everything for that grid point in a pickle file
    '''
    extra_pickle_extension = '/home/gastonlagaffe/../../vol_c'
    pickle.dump( psf_stats, open( extra_pickle_extension + '/pickle_files/chrom_' + chromaticString + '_avgwavel_' + wavelString +
                                  '_opd_' + opdString + '_tip_' + tipString + '_tilt_' + tiltString +
                                  '_transl_' + translString + '_PS_' + PSstring + extraString + '.pkl', "wb" ) )
    '''

    # save a simple FITS file of the PSF, along with aberration info in the header
    #fits_extension = "~/../../Volumes/seagate_external_drive/lbti_data_reduction/190526_fake_fizeau/"
    fits_extension = "./junk/"
    
    img_size = np.shape(psf_stats['PSF_image'])[0]
    half_cutout = 1024 # this sets the size of the saved image
    simple_fits = psf_stats['PSF_image'][int(0.5*img_size)-half_cutout:int(0.5*img_size)+half_cutout,
                                         int(0.5*img_size)-half_cutout:int(0.5*img_size)+half_cutout]
    simple_fits = np.float32(simple_fits)

    # also include FFT amplitude and phase
    amp_fft = psf_stats['FTamp_image'][int(0.5*img_size)-half_cutout:int(0.5*img_size)+half_cutout,
                                         int(0.5*img_size)-half_cutout:int(0.5*img_size)+half_cutout]
    phase_fft = psf_stats['FTphase_image'][int(0.5*img_size)-half_cutout:int(0.5*img_size)+half_cutout,
                                         int(0.5*img_size)-half_cutout:int(0.5*img_size)+half_cutout]

    cube_to_write = np.ones((3,np.shape(simple_fits)[0],np.shape(simple_fits)[1]))
    cube_to_write[0,:,:] = simple_fits
    cube_to_write[1,:,:] = amp_fft
    cube_to_write[2,:,:] = phase_fft

    #hdu = fits.PrimaryHDU(simple_fits)

    hdr = fits.Header()
    hdr["WAVEL_UM"] = np.float(avg_wavel*1e6) # wavelength (um)
    hdr["MONOCHR"] = monochromatic # monochromatic? (true/false)
    hdr["OPD_UM"] = np.float(opd*1e6) # OPD (um)
    hdr["TIPY_MAS"] = np.float(diff_tip*1e3) # tip (mas)
    hdr["TILTXMAS"] = np.float(diff_tilt*1e3) # tilt (mas)


    # the following FITS command was written with random walk PSFs in mind
    #fits.writeto(fits_extension + "psf_trial1_" + str("{:0>8d}".format(int(index_num))) + ".fits", data=simple_fits, header=hdr, overwrite=True)
    # ... the alternative is to just give it an ersatz name
    fits.writeto("test.fits", data=cube_to_write, header=hdr, overwrite=True)

    
    #plt.plot()
    #plt.savefig("junk.png")

    elapsed_time = time.time() - start_time
    print('Elapsed time for this PSF: '+str(elapsed_time))


####################

def main():

    ######################
    ## USER INPUTS

    mode = "permutations"
    #mode = "random_walk"
    
    opd_start = 0.0e-6
    #opd_stop = 0.0e-6
    opd_stop = 1.0e-6 # inclusive
    opd_increment = 3.0e-6 # change in OPD at each step; in m

    tilt_start = 0.0
    tilt_stop = 1.0 # asec 
    tilt_increment = 2.0

    tip_start = 0.0
    tip_stop = 1.0 # asec 
    tip_increment = 2.0 

    transl_start = 0.0 # position at which to start
    transl_stop = 1.0
    transl_increment = 3.0

    extra_pickle_title_string = 'junk_psf_stats'

    ## END USER INPUTS
    ######################

    start_time_very_beginning = time.time()

    if (mode == "permutations"):
        '''
        Generate all possible permutations of the input min, max, and step
        sizes of parameters
        '''
        
        # generate the grid arrays we want
        num_opd = np.divide(opd_stop-opd_start, opd_increment) + 1 # number of samples
        opdArray = np.linspace(opd_start, opd_stop, num=int(num_opd), endpoint=True)
        num_tip = np.divide(tip_stop-tip_start, tip_increment) + 1
        tipArray = np.linspace(tip_start, tip_stop, num=int(num_tip), endpoint=True)
        num_tilt = np.divide(tilt_stop-tilt_start, tilt_increment) + 1
        tiltArray = np.linspace(tilt_start, tilt_stop, num=int(num_tilt), endpoint=True)
        num_transl = np.divide(transl_stop-transl_start, transl_increment) + 1
        translArray = np.linspace(transl_start, transl_stop, num=int(num_transl), endpoint=True)

        # put everything into a list of lists of combinations
        combineArray = [opdArray,tipArray,tiltArray,translArray]
        permutationsArray = list(itertools.product(*combineArray))

        # this is just a kludge needed to make 'permutations' mode work, after implementation
        # of 'random_walk'
        index_array = [[1]]

        params_array = np.concatenate((permutationsArray,index_array), axis=1)

    elif (mode == "random_walk"):
        '''
        Generate parameters as if they were making a random walk
        '''

        opd_scale = 50e-6 # m
        tip_y_scale = 0.2 # asec
        tilt_x_scale = 0.2 # asec

        # change these as necessary
        N = 10000
        opdArray = walker(opd_scale, N)
        #opdArray = np.zeros(N)
        tipArray = walker(tip_y_scale, N)
        #tipArray = np.zeros(N)
        tiltArray = walker(tilt_x_scale, N)
        #tiltArray = np.zeros(N)
        translArray = np.zeros(N+1)

        index_array = np.arange(0,N+1)

        print("Random walk, OPD array head:")
        print(opdArray[0:10])
        print("Random walk, tip array head:")
        print(tipArray[0:10])
        print("Random walk, tilt array head:")
        print(tiltArray[0:10])

        plt.plot(np.multiply(1e6,opdArray))
        plt.title("OPD (um)")
        plt.savefig("opd_ref.pdf")
        plt.clf()

        plt.plot(tipArray)
        plt.title("Tip_y (asec)")
        plt.savefig("tip_ref.pdf")
        plt.clf()

        plt.plot(tiltArray)
        plt.title("Tilt_x (asec)")
        plt.savefig("tilt_ref.pdf")
        plt.clf()

        # put everything into a list of lists of combinations
        combineArray = np.stack((opdArray,tipArray,tiltArray,translArray,index_array), axis=0).T
        params_array = combineArray
        

    # read in masks (0=masked; 1=good) and change convention (True=masked; False=good)
    global mask_all
    global mask_center
    global mask_right
    mask_all0 = fits.open('masks/mask_all_power.fits')
    mask_all = np.abs(np.subtract(mask_all0[0].data,1)).astype(bool)
    mask_center0 = fits.open('masks/mask_center_power.fits')
    mask_center = np.abs(np.subtract(mask_center0[0].data,1))
    mask_right0 = fits.open('masks/mask_right_power.fits')
    mask_right = np.abs(np.subtract(mask_right0[0].data,1)).astype(bool)
    
    ncpu = mp.cpu_count()

    # print job info to screen
    print('----------')
    print('Number cores: ')
    print(ncpu)
    print('----------')
    print('Commencing grid job with ')
    print('----------')
    print('OPD of ')
    print(opdArray)
    print('----------')
    print('Tip of')
    print(tipArray)
    print('----------')
    print('Tilt of')
    print(tiltArray)
    print('----------')
    print('Translation of')
    print(translArray)
    print('----------')
        
    pool = Pool(ncpu) # create pool object
    list_dicts = pool.map(psf_sim, params_array)

    elapsed_time_since_beginning = time.time() - start_time_very_beginning
    print('Total elapsed time: '+str(elapsed_time_since_beginning))

##################

# do it!
if __name__ == '__main__':
    main()
