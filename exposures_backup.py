from lmircam_tools import *
from lmircam_tools.print_tools import info, request

def get_lmircam_frames(dit, coadds, nseqs, use_bg=0, save_data=True):           #  Standard function to execute LMIRCam exosures.
    pi.setINDI("LMIRCAM.stop.value=Off")                                        #  Make sure the camera is not running.
    
    if save_data:                                                               #  If data saving is enabled.
        pi.setINDI("LMIRCAM.enable_save.value=On")                              #  Turn on data saving.
    else:                                                                       #  Else.
        pi.setINDI("LMIRCAM.enable_save.value=Off")                             #  Turn off data saving.
    
    #  Now run the camera with the parameters provided:
    pi.setINDI("LMIRCAM.acquire.int_time=%f;num_coadds=%i;num_seqs=%i;enable_bg=%i;is_bg=0;is_cont=0"%(dit, coadds, nseqs, use_bg), timeout = 4.0 * dit * nseqs)
    # Really just guessing on that timeout in the above line.  The way I understand this, this timeout should not matter, and most of the time it doesn't.  But then in some cases I get a timeout error if I don't set this.

#    try:  #  See if it works.
#        pi.setINDI("LMIRCAM.acquire.int_time=%f;num_coadds=%i;num_seqs=%i;enable_bg=%i;is_bg=0;is_cont=0"%(dit, coadds, nseqs, use_bg), timeout = 2.0 * dit * nsequ)
#    except:  #  What if not?  Probably no background available of available backgroun invalid?  Try without background.
#        info('Could not start exposure, trying without background.')
#        try:  #  See if this works now.
#            pi.setINDI("LMIRCAM.acquire.int_time=%f;num_coadds=%i;num_seqs=%i;enable_bg=0;is_bg=0;is_cont=0"%(dit, coadds, nseqs, use_bg), timeout = 2.0 * dit * nsequ)
#        except:  #  If it still doesn't work?  Some deeper problem that needs attention.  Tell operator and give up.
#            info('Still unable to start exposure.  Giving up.')
#            request('Please fix problem with exposure and restart.')
#            quit()


def get_lmircam_bg(dit, coadds, save_data=True):                                #  Standard function to take an LMIRCam background (saving or not saving the frame).
    pi.setINDI("LMIRCAM.stop.value=Off")                                        #  Make sure the camera is not running.
    
    if save_data:                                                               #  If data saving is enabled.
        pi.setINDI("LMIRCAM.enable_save.value=On")                              #  Turn on data saving.
    else:                                                                       #  Else.
        pi.setINDI("LMIRCAM.enable_save.value=Off")                             #  Turn off data saving.

    #  Now run the camera with the parameters provided:
    pi.setINDI("LMIRCAM.acquire.int_time=%f;num_coadds=%i;num_seqs=1;enable_bg=0;is_bg=1;is_cont=0"%(dit, coadds), timeout = 4.0 * dit)
