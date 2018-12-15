import numpy as np
from numpy import ma
from pyindi import *

pi = PyINDI(verbose=False)

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

def put_in_grism():
    '''
    Inserts the LMIR grism
    '''
    pi.setINDI("Lmir.lmir_FW3.command", "Lgrism6AR", timeout=45, wait=True)
    pi.setINDI("Lmir.lmir_FW25.command", "Lspec2.8-4.0", timeout=45, wait=True) # blocks some extraneous light

    
def remove_grism():
    '''
    Inserts the LMIR grism
    '''
    pi.setINDI("Lmir.lmir_FW3.command", "## WHATEVER FILTER GOES HERE ##", timeout=45, wait=True)
    pi.setINDI("Lmir.lmir_FW25.command", "## WHATEVER FILTER GOES HERE ##", timeout=45, wait=True) # blocks some extraneous light
    
        
'''
def check_ao_pc_loops():
    # check that the AO and Phasecam loops are closed

    ## ## INDI COMMAND TO SEE LEFT AO STATUS, RIGHT AO STATUS


class check_pc_loop(check_ao_loops):
    # check that Phasecam loop is closed (inherits check of AO loops)
    
    ## ## INDI COMMAND TO SEE PHASECAM STATUS

    ## ## IF AO AND PC ARE ALL CLOSED, RETURN 2
    ## ## IF AO ARE CLOSED BUT PC IS NOT, RETURN 1

'''
