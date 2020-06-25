# This tests Fizeau aberration retrieval algorithms locally

import pdb
from lmircam_tools import *
from lmircam_tools.overlap_psfs import overlap_psfs
from lmircam_tools.dial_opd import find_optimal_opd_fizeau_grism, live_opd_correction_fizeau_grism, implement_optimal_opd
from lmircam_tools.change_tt import print_write_fft_info, get_apply_pc_setpts, compare_setpts

# print fft info, see how it compares with the set thresholds
# modes: "fake_fits" / "az_source" / "science"
num_psfs, fftimg_shape = print_write_fft_info(integ_time, sci_wavel = wavel_lambda, mode = mode_choice)

# calculate and apply Phasecam setpoints; write them to a text file to check the correction
setpoints_pickle_pre = "setpoints_pickle_pre.pkl"

# function to read files from disk, analyze PSFs, and
# 1. calculate TT and OPD
# 2. write results to file, along with injected values from FITS headers
retrieve_and_write()

get_apply_pc_setpts(integ_time,
					num_psfs,
					fftimg_shape,
					sci_wavel = wavel_lambda,
					mode = mode_choice,
					pickle_name = setpoints_pickle_pre,
					apply = True)
