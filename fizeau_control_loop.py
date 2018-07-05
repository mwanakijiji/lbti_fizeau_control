#!/usr/bin/python

from lmircam_tools.overlap_psfs import overlap_airy_psfs
from lmircam_tools.


############## BEGIN GROSS OVERLAP OF AIRY PSFS

psf_loc_setpoint = [1220,800]  # pixel location for PSFs to be at
overlap_airy_psfs(psf_loc_setpoint)

############## END GROSS OVERLAP OF AIRY PSFS



############## BEGIN PUT IN GRISM AND REFINE GRISM-PSF OVERLAP

overlap_grism_psfs(psf_loc_setpoint)

############## END PUT IN GRISM AND REFINE GRISM-PSF OVERLAP


############## BEGIN DIAL OPD WITH HPC AND FIND CENTER OF COHERENCE ENVELOPE, THEN REMOVE GRISM

dial_opd_fizeau_grism(psf_loc_setpoint) # might also use argument of the re-established Fizeau/grism PSF instead of the coordinate where it's supposed to be

############## END DIAL OPD WITH HPC AND FIND CENTER OF COHERENCE ENVELOPE, THEN REMOVE GRISM=


############## BEGIN HOLD CENTER OF SCIENCE COHERENCE ENVELOPE WITH HIGH-CONTRAST FRINGES

# maybe allow HPC to correct pathlength using science readouts, even though Phasecam not closed yet? or does this not have to be automatic, and we can just correct the Fizeau/Airy PSFs once Phasecam is closed?

############## END HOLD CENTER OF SCIENCE COHERENCE ENVELOPE WITH HIGH-CONTRAST FRINGES


############## TRANSLATE NIL + CLOSE PHASECAM LOOP HERE?
# might be manual step


############## BEGIN OPTIMIZE SCIENCE PSF BY FINDING OPD AND TT SETPOINTS ITERATIVELY

dial_opd_fizeau_airy(psf_location) 
## adjust TT to optimize PSF; maybe iterate with OPD?
## note OPD movements cannot be more than 5 um with Phasecam closed

############## END OPTIMIZE SCIENCE PSF BY FINDING OPD AND TT SETPOINTS ITERATIVELY


############## BEGIN SET FLAG IN SCIENCE IMAGE HEADERS TO EFFECT FIZEAU IS CONTROLLED


############## END SET FLAG IN SCIENCE IMAGE HEADERS TO EFFECT FIZEAU IS CONTROLLED



############## KEEP TRACK OF FRINGES IN SCIENCE READOUTS



############## REDO EVERYTHING ONCE NOD HAPPENS
