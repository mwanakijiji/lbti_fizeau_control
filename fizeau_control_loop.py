#!/usr/bin/python

import pdb
from lmircam_tools import *
from lmircam_tools.overlap_psfs import overlap_psfs
from lmircam_tools.dial_opd import find_optimal_opd_fizeau_grism, live_opd_correction_fizeau_grism, implement_optimal_opd
from lmircam_tools.change_tt import print_write_fft_info, get_apply_pc_setpts, compare_setpts

#################################################################################
###############         ALIGNMENT         #######################################
#################################################################################

############## GROSS OVERLAP OF NON-FIZEAU AIRY PSFS

# ROI size no larger than 512x512, but large enough in y to cover 200-pixel-tall
#   grism PSF + possible movement; say, ~400x200
# insert option to only horizontally overlap grisms

# default integration time is set in the __init__ file (vestigial now --2019Feb10)

print("----------------------------------------------------------------------------------")
raw_input("Place an ROI, no larger than 512x512, over the Phasecam sweet spot for LMIRcam, and " + \
	"move telescopes\n to put Airy PSFs somethere inside the ROI. \nROI size needs to be "+\
	"large enough to include the grism PSF and possible movement\nPress [Enter] when done.")

# Operating mode
#    [MODE]           [SOURCE OF DATA]               [DATA CONFIGURATION]  [CMDS TO MIRRORS AND WHEELS]  [CMDS TO TELESCOPE]
#    "fake_fits":     fake/old FITS files            written file read in  N                             N
#    "az_source":     piggy-back off frames user     INDI                  Y                             N
#                       is taking with LMIR
#    "science":       (same as "az_source")          INDI                  Y                             Y
print("----------------------------------------------------------------------------------")
mode_choice = "az_source"
print("This optimization code is running in mode " + mode_choice)
print("Stop continuous aquisition of the camera for alignment and set the integration time in the init file.")
print("----------------------------------------------------------------------------------")

print("Science wavelength:")
print(wavel_lambda)

pdb.set_trace()
# modes: "az_source" / "science"
overlap_psfs(integ_time, fiz_lmir_sweet_spot, mode = mode_choice, psf_type = "airy") # filter-agnostic


############## PUT IN GRISM AND REFINE GRISM-PSF OVERLAP

# modes: "az_source" / "science"
put_in_grism(mode = mode_choice)

# modes: "az_source" / "science"
overlap_psfs(integ_time, fiz_lmir_sweet_spot, mode = mode_choice, psf_type = "grism")


############## IN GRISM MODE, DIAL OPD WITH HPC AND FIND CENTER OF COHERENCE ENVELOPE, THEN REMOVE GRISM

# modes: "fake_fits" / "az_source" / "science"
# this should be run WHILE taking science data
live_opd_correction_fizeau_grism(integ_time, mode = mode_choice)

# modes: "fake_fits" / "az_source" / "science"
find_optimal_opd_fizeau_grism(integ_time, mode = mode_choice)
implement_optimal_opd(mode = mode_choice)
print("----------------------------------------------------------------------------------")
raw_input("Science detector alignment done. Now align Phasecam and close the phase loop")


#################################################################################
###############         SCIENCE           #######################################
#################################################################################


# print fft info, see how it compares with the set thresholds
# modes: "fake_fits" / "az_source" / "science"
num_psfs, fftimg_shape = print_write_fft_info(integ_time, sci_wavel = wavel_lambda, mode = mode_choice)

# calculate and apply Phasecam setpoints; write them to a pickle file to check the correction
setpoints_pickle_pre = "setpoints_pickle_pre.pkl"
# modes: "fake_fits" / "az_source" / "science"
get_apply_pc_setpts(integ_time, num_psfs, fftimg_shape, sci_wavel = wavel_lambda, mode = mode_choice, pickle_name = setpoints_pickle_pre, apply = True)

# print/calculate FFT info after the change so as to check setpoints for sign (code cant tell which PSF is SX and DX)
# modes: "fake_fits" / "az_source" / "science"
num_psfs, fftimg_shape = print_write_fft_info(integ_time, sci_wavel = wavel_lambda, mode = mode_choice, checker=True) # note checker=True, because we're checking the correction

# recalculate setpoints, but dont apply a new correction just yet
setpoints_pickle_post = "setpoints_pickle_post.pkl"
# modes: "fake_fits" / "az_source" / "science"
get_apply_pc_setpts(integ_time, num_psfs, fftimg_shape, sci_wavel = wavel_lambda, mode = mode_choice, pickle_name = setpoints_pickle_post, apply = False)

# compare setpoints; if any one aspect (TT, PL) is worse, flip it back the other way, 2x
# modes: "fake_fits" / "az_source" / "science"
compare_setpts(setpoints_pickle_pre, setpoints_pickle_post, mode = mode_choice)
