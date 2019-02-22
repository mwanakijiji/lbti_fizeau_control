from pyindi import *
from lmircam_tools.utils import nod
import numpy as np

integ_time = 0.1

pi = PyINDI(host="localhost", verbose=False)

pi.setINDI("fizeau.enable_run.value=On")

#pi.setINDI("LMIRCAM_save.enable_save.value=On")

f = pi.getFITS("fizeau.roi_image.file", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=5)

# write roi_image to file
image = f[0].data
hdu = pyfits.PrimaryHDU(image)
hdulist = pyfits.HDUList([hdu])
hdu.writeto("roi_image.fits", clobber=True)

''' 
f = pi.getFITS("fizeau.mag_image.file", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=5)

# write roi_image to file
image = f[0].data
hdu = pyfits.PrimaryHDU(image)
hdulist = pyfits.HDUList([hdu])
hdu.writeto("mag_image.fits", clobber=True)

f = pi.getFITS("fizeau.phase_image.file", "LMIRCAM.acquire.enable_bg=1;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % integ_time, timeout=5)

# write roi_image to file
image = f[0].data
hdu = pyfits.PrimaryHDU(image)
hdulist = pyfits.HDUList([hdu])
hdu.writeto("phase_image.fits", clobber=True)
'''
