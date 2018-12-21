#!/usr/bin/python

import pdb
from lmircam_tools import *
from lmircam_tools.overlap_psfs import overlap_psfs
from lmircam_tools.dial_opd import optimize_opd_fizeau_grism, optimize_opd_fizeau_airy
from lmircam_tools.change_tt import get_pc_setpts

############## GROSS OVERLAP OF NON-FIZEAU AIRY PSFS

# ROI size no larger than 512x512, but large enough in y to cover 200-pixel-tall
# grism PSF + possible movement; ~400x200
# insert option to only horizontally overlap grisms

print("----------------------------------------------------------------------------------")
raw_input("Place an ROI, no larger than 512x512, over the Phasecam sweet spot for LMIRcam, and move telescopes\n to put Airy PSFs somethere inside the ROI. \nROI size needs to be "+\
	"large enough to include the grism PSF and possible movement, e.g., ~400x200\n (HOSTS sweet spot was [880,800] "+\
	"on full detector in May 2018.)\nPress [Enter] when done.")

# Phasecam sweet spot on detector, in ROI coordinates (y, x)
fiz_lmir_sweet_spot = [200,100] 

# Is this script being run as a test, or are we doing on-sky science?
# Options: "fake_fits", "science"
print("----------------------------------------------------------------------------------")
mode_choice = "fake_fits"
print("This optimization code is running in mode " + mode_choice)

pdf.set_trace()
overlap_psfs(fiz_lmir_sweet_spot, mode = mode_choice, psf_type = "airy") # filter-agnostic

## ## see old sweet spots (can also locate them on NOMIC, and then see where they are on LMIR)
## ## see nomic nulling to see how nod with wheel is done

############## PUT IN GRISM AND REFINE GRISM-PSF OVERLAP

put_in_grism()
overlap_psfs(fiz_lmir_sweet_spot, mode = mode_choice, psf_type = "grism")


############## IN GRISM MODE, DIAL OPD WITH HPC AND FIND CENTER OF COHERENCE ENVELOPE, THEN REMOVE GRISM

## ## remove psf_loc_setpoint from func; obsolete optim airy fcn
## ## add length of scan (in total OPD?)
optimize_opd_fizeau_grism(mode = mode_choice) # might also use argument of the re-established Fizeau/grism PSF instead of the coordinate where it's supposed to be
remove_grism()
print("----------------------------------------------------------------------------------")
raw_input("Now align Phasecam and close the phase loop")


############## HOLD CENTER OF SCIENCE COHERENCE ENVELOPE WITH HIGH-CONTRAST FRINGES

# 2ND PRIORITY
# maybe allow HPC to correct pathlength using science readouts, even though Phasecam not closed yet? or does this not have to be automatic, and we can just correct the Fizeau/Airy PSFs once Phasecam is closed?


############## TRANSLATE NIL + CLOSE PHASECAM LOOP HERE?

# 2ND PRIORITY
# might be manual step


############## OPTIMIZE SCIENCE PSF BY FINDING OPD AND TT SETPOINTS ITERATIVELY

## ## optimize_opd_fizeau_airy(psf_location) # 2ND PRIORITY
get_pc_setpts(log_name = "setpt_log.csv", mode = mode_choice) 
## adjust TT to optimize PSF; maybe iterate with OPD?
## note OPD movements cannot be more than 5 um with Phasecam closed


############## ANY RUNNING STATS I WANT TO KEEP OF THE SCIENCE PSFS?



############## REDO EVERYTHING ONCE NOD HAPPENS
