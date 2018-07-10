import time
## ## from pyindi import *
import numpy as np
from numpy import ma

## ## pi = PyINDI(verbose=False)



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
