#!/usr/bin/python

from pyindi import *
pi = PyINDI()

print "getting display image"
f = pi.getFITS("LMIRCAM.DisplayImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)
#print f[0].header["NAXIS1"], "x", f[0].header["NAXIS2"]
'''
print "turning on fizRun";
pi.setINDI ("LMIRCAM.fizRun.value=On");

print "getting PSF";
f = pi.getFITS("LMIRCAM.fizPSFImage.File", "LMIRCAM.acquire.enable_bg=0;int_time=%i;is_bg=0;is_cont=0;num_coadds=1;num_seqs=1" % 100, timeout=60)
print f[0].header["NAXIS1"], "x", f[0].header["NAXIS2"]

'''
