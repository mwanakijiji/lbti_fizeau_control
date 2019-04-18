from lmircam_tools import *
from lmircam_tools.print_tools import info, request
from time import sleep
inner_loop_max = 100
sleep_time = 10.0  # sleep time between N_sequences exposures


def get_lmircam_frames(dit, coadds, nseqs, use_bg=0, inner_loop_max=inner_loop_max, save_data=True):           #  Standard function to execute LMIRCam exosures.
    pi.setINDI("lmircam.stop.value=Off")                                        #  Make sure the camera is not running.
    
    N_outer_loops = int(nseqs) / int(inner_loop_max)
    remainder = nseqs % inner_loop_max
    
    if use_bg==1:
      pi.setINDI("lmircam_display.enable_bg.value=On")
    else:
      pi.setINDI("lmircam_display.enable_bg.value=Off")
    
    if save_data:                                                               #  If data saving is enabled.
        pi.setINDI("lmircam_save.enable_save.value=On")                              #  Turn on data saving.
    else:                                                                       #  Else.
        pi.setINDI("lmircam_save.enable_save.value=Off")                             #  Turn off data saving.
    
    #  Now run the camera with the parameters provided:
    for ii in xrange(N_outer_loops):
        pi.setINDI("lmircam.acquire.int_time=%f;is_cont=0;num_coadds=%i;num_seqs=%i"%(dit, coadds, inner_loop_max), timeout = 100000000.0)
        sleep(sleep_time)
    if remainder > 0:
        pi.setINDI("lmircam.acquire.int_time=%f;is_cont=0;num_coadds=%i;num_seqs=%i"%(dit, coadds, remainder), timeout = 100000000.0)

# OUTDATED, don't use, won't work.  Only left for future reference.
#def get_lmircam_bg(dit, coadds, save_data=True):                                #  Standard function to take an LMIRCam background (saving or not saving the frame).
#    pi.setINDI("LMIRCAM.stop.value=Off")                                        #  Make sure the camera is not running.
#    
#    if save_data:                                                               #  If data saving is enabled.
#        pi.setINDI("LMIRCAM_save.enable_save.value=On")                              #  Turn on data saving.
#    else:                                                                       #  Else.
#        pi.setINDI("LMIRCAM_save.enable_save.value=Off")                             #  Turn off data saving.
#
#    #  Now run the camera with the parameters provided:
#    pi.setINDI("LMIRCAM.acquire.int_time=%f;num_coadds=%i;num_seqs=1;enable_bg=0;is_bg=1;is_cont=0"%(dit, coadds), timeout = 100000000.0)
