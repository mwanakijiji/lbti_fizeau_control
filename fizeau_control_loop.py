#!/usr/bin/python

import pdb
from lmircam_tools import *
from lmircam_tools.overlap_psfs import overlap_psfs
from lmircam_tools.dial_opd import find_optimal_opd_fizeau_grism, live_opd_correction_fizeau_grism, implement_optimal_opd
from lmircam_tools.change_tt import print_write_fft_info, get_apply_pc_setpts

############## GROSS OVERLAP OF NON-FIZEAU AIRY PSFS

# ROI size no larger than 512x512, but large enough in y to cover 200-pixel-tall
# grism PSF + possible movement; ~400x200
# insert option to only horizontally overlap grisms

# default integration time is set in the __init__ file

print("----------------------------------------------------------------------------------")
raw_input("Place an ROI, no larger than 512x512, over the Phasecam sweet spot for LMIRcam, and move telescopes\n to put Airy PSFs somethere inside the ROI. \nROI size needs to be "+\
	"large enough to include the grism PSF and possible movement\nPress [Enter] when done.")

# Phasecam sweet spot on detector, in ROI coordinates (y, x)
fiz_lmir_sweet_spot = [200,100] 

# Is this script being run as a test, or are we doing on-sky science?
# Options: "total_passive", "fake_fits", "artif_source", "science"
#    "total_passive": read in fake FITS files, and no getFITS or setINDI commands are sent (but getINDI are)
#    "fake_fits":     read in fake FITS files, but also send INDI commands to cameras and mirrors (but not the telescope)
#    "artif_source":  use detector images involving artificial sources in closed-dome, and send commands to UBC mirrors (but not the telescope)
#    "science":       send commands to cameras, mirrors, and telescope like we're on-sky
print("----------------------------------------------------------------------------------")
mode_choice = "fake_fits"
print("This optimization code is running in mode " + mode_choice)
print("Stop continuous aquisition of the camera.")
print("----------------------------------------------------------------------------------")

pdb.set_trace()
overlap_psfs(integ_time, fiz_lmir_sweet_spot, mode = mode_choice, psf_type = "airy") # filter-agnostic

## ## see old sweet spots (can also locate them on NOMIC, and then see where they are on LMIR)
## ## see nomic nulling to see how nod with wheel is done

############## PUT IN GRISM AND REFINE GRISM-PSF OVERLAP

put_in_grism()
overlap_psfs(integ_time, fiz_lmir_sweet_spot, mode = mode_choice, psf_type = "grism")


############## IN GRISM MODE, DIAL OPD WITH HPC AND FIND CENTER OF COHERENCE ENVELOPE, THEN REMOVE GRISM

## ## remove psf_loc_setpoint from func; obsolete optim airy fcn
## ## add length of scan (in total OPD?)
## ## sometimes new HPC movement causes grisms to separate; may need to re-overlap them each time
## ## insert all hpc, fpc piston and TT statuses into headers

# the following function is for a one-off correction once you have a grism PSF with visible fringes at an angle
# once its seen to work well, then it should be applied periodically while data is being taken (like once the angle is see to be >5 degrees, or something like that)
## ## ONLY WORKS IN GRISM MODE, THOUGH
live_opd_correction_fizeau_grism(integ_time, mode = mode_choice)

# the following is for doing a big, automated scan... but it may be more efficient to do the scan manually
find_optimal_opd_fizeau_grism(integ_time, mode = mode_choice) # might also use argument of the re-established Fizeau/grism PSF instead of the coordinate where it's supposed to be
implement_optimal_opd(mode = mode_choice)
print("----------------------------------------------------------------------------------")
raw_input("Now align Phasecam and close the phase loop")

## ## ... AND REMOVE THE GRISM TOO

############## HOLD CENTER OF SCIENCE COHERENCE ENVELOPE WITH HIGH-CONTRAST FRINGES

# 2ND PRIORITY
# maybe allow HPC to correct pathlength using science readouts, even though Phasecam not closed yet? or does this not have to be automatic, and we can just correct the Fizeau/Airy PSFs once Phasecam is closed?


############## TRANSLATE NIL + CLOSE PHASECAM LOOP HERE?

# 2ND PRIORITY
# might be manual step


############## OPTIMIZE SCIENCE PSF BY FINDING OPD AND TT SETPOINTS ITERATIVELY

# print fft info, see how it compares with the set thresholds
print_write_fft_info(integ_time, mode = mode_choice, fft_pickle_write_name = "fft_info.pkl")

# calculate and apply Phasecam setpoints
## ## DONT RUN THE FOLLOWING FUNCTION YET; JUST USE BITS INSIDE OF IT TO TEST THE EFFECT ON THE PSF
get_apply_pc_setpts(integ_time, mode = mode_choice, fft_pickle_read_name = "fft_info.pkl")

## adjust TT to optimize PSF; maybe iterate with OPD?
## note OPD movements cannot be more than 5 um with Phasecam closed


############## ANY RUNNING STATS I WANT TO KEEP OF THE SCIENCE PSFS?



############## REDO EVERYTHING ONCE NOD HAPPENS
