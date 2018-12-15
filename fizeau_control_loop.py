#!/usr/bin/python

import pdb
from lmircam_tools import *
from lmircam_tools.overlap_psfs import overlap_airy_psfs, overlap_grism_psfs
from lmircam_tools.dial_opd import optimize_opd_fizeau_grism, optimize_opd_fizeau_airy
from lmircam_tools.change_tt import optimize_tt_fizeau_airy


############## BEGIN GROSS OVERLAP OF NON-FIZEAU AIRY PSFS

psf_loc_setpoint = [150,40]  # pixel location for PSFs to be at (y,x), zeroed to the ROI region
overlap_airy_psfs(psf_loc_setpoint) # filter-agnostic

############## END GROSS OVERLAP OF AIRY PSFS



############## BEGIN PUT IN GRISM AND REFINE GRISM-PSF OVERLAP

#put_in_grism()
#overlap_grism_psfs(psf_loc_setpoint)

############## END PUT IN GRISM AND REFINE GRISM-PSF OVERLAP


############## BEGIN DIAL OPD WITH HPC AND FIND CENTER OF COHERENCE ENVELOPE, THEN REMOVE GRISM

#optimize_opd_fizeau_grism(psf_loc_setpoint) # might also use argument of the re-established Fizeau/grism PSF instead of the coordinate where it's supposed to be
#remove_grism()

############## END DIAL OPD WITH HPC AND FIND CENTER OF COHERENCE ENVELOPE, THEN REMOVE GRISM


############## BEGIN HOLD CENTER OF SCIENCE COHERENCE ENVELOPE WITH HIGH-CONTRAST FRINGES

# 2ND PRIORITY
# maybe allow HPC to correct pathlength using science readouts, even though Phasecam not closed yet? or does this not have to be automatic, and we can just correct the Fizeau/Airy PSFs once Phasecam is closed?

############## END HOLD CENTER OF SCIENCE COHERENCE ENVELOPE WITH HIGH-CONTRAST FRINGES


############## TRANSLATE NIL + CLOSE PHASECAM LOOP HERE?

# 2ND PRIORITY
# might be manual step


############## BEGIN OPTIMIZE SCIENCE PSF BY FINDING OPD AND TT SETPOINTS ITERATIVELY

## ## optimize_opd_fizeau_airy(psf_location) # 2ND PRIORITY
#optimize_tt_fizeau_airy(psf_location) 
## adjust TT to optimize PSF; maybe iterate with OPD?
## note OPD movements cannot be more than 5 um with Phasecam closed

############## END OPTIMIZE SCIENCE PSF BY FINDING OPD AND TT SETPOINTS ITERATIVELY


############## ANY RUNNING STATS I WANT TO KEEP OF THE SCIENCE PSFS?



############## REDO EVERYTHING ONCE NOD HAPPENS
