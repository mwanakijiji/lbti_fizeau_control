import sys
sys.path.append('/home/observer/Scripts/observing_scripts/tests/')
import time
from pyindi import * 
import numpy as np
from lmircam_tools.utils import wait4AORunning, setFLAG
from lmircam_tools.exposures import get_lmircam_frames

pi = PyINDI(verbose = False)

def ten(coord):
    '''change a sexageismal to a decimal
    INPUTS:
    coord a python list [hh/dd, mm, ss]'''
    coord=np.array(coord)
    num=len(coord)
    fac=np.array([1.,60.0,3600.0])
    if (coord < 0).sum() > 0:
        sign=-1
    else:
        sign=1
    return sign*(np.abs(coord)/fac[:num]).sum()

def makeOffsets(RA,Dec,sep,PA):
    '''Take the RA and Dec of the primary, and the 
    separation and PA of the secondary, and return
    the appropriate nod vectors for 3-point dither
    pattern.
    PA East of North. North up East left.
    INPUTS:
    RA = [hh, mm, ss.dd]
    Dec = [dd, mm, ss.dd]
    sep = arcseconds
    PA = degrees'''
    secDecOff = -1*sep*np.cos(PA*np.pi/180.)#arcseconds
    secRAOff = -1*sep*np.sin(PA*np.pi/180.)#arcseconds NO cos(dec)!!!

    skyDecOff = 2*secDecOff
    skyRAOff = 2*secRAOff

    priDecOff = 0
    priRAOff = 0

    return priRAOff, priDecOff, secRAOff, secDecOff, skyRAOff, skyDecOff
    
def nod_sneakSomeDarks(RAOff,DecOff,darkIntParams,side,coords="RADEC",absrel='ABS'):
    '''While moving
    spin lmir_fw2 by half a position to act as a blank and grab some
    dark frames. Then move lmir_fw2 back to where it was.'''
####  Here are some lines to fix DX issue.
    if side == 'right' or side == 'both':    # This is only relevant if DX is involved.
      low_gain_file = 'gain0.1_10modi.fits'  # File name for low gain file.
      settling_time = 3                      # Time the AO should be closed in low gains before setting high gain again.
      full_gain_path = pi.getINDI('indi_dx_adsec_msgd.ADSEC___R___G_GAIN_A.value')  # Path of the full gain file.
      low_gain_path = os.path.join(os.path.dirname(full_gain_path), low_gain_file)  # Path of the low gain file (same directory, but above file name).
      print 'Setting 10-mode gain.'          # Some status info.
      pi.setINDI('ao_indi_dx_ice.load_gain_tag_a.tag', low_gain_path, wait=True)    # Set low gains before nodding.
####  Those were some lines to fix DX issue.
    pi.setINDI("LBTO.OffsetPointing.CoordSys", coords,
           "LBTO.OffsetPointing.OffsetX", DecOff, 
           "LBTO.OffsetPointing.OffsetY", RAOff, 
           "LBTO.OffsetPointing.Side", side, 
           "LBTO.OffsetPointing.Type", absrel, 
           timeout=400,
           wait=False) 
	###where's FW2
    lightIntParams = pi.getINDI("LMIRCAM.acquire")
    currentFlag = pi.getINDI("LMIRCAM.EditFITS")
    setFLAG('DRK')
    fw2_pos = pi.getINDI("Lmir.lmir_FW2_status.PosNum", wait=True)
    pi.setINDI("Lmir.lmir_FW2.command", int(fw2_pos)+25000, timeout=20, wait=True)
    get_lmircam_frames(*darkIntParams,use_bg=True)
    pi.setINDI("Lmir.lmir_FW2.command", int(fw2_pos), timeout=20, wait=True)
    pi.setINDI(currentFlag)
    wait4AORunning(side)
####  Here are some more lines to fix DX issue.
    if side == 'right' or side == 'both':  # This is only relevant if DX is involved.
      print 'Waiting a few seconds...'     # Some status info.
      sleep(settling_time)                 # Sleep to allow AO to settle in low gain.
      pi.setINDI('ao_indi_dx_ice.load_gain_tag_a.tag', full_gain_path, wait=True)  # Set gains back to original, high values.
      print 'Ok to continue.'              # All done.
####  Those were some more lines to fix DX issue.

def nod_sneakSomeDarks_normal(RAOff,DecOff,darkIntParams,side,coords="RADEC",absrel='ABS'):
    '''While moving
    spin lmir_fw2 by half a position to act as a blank and grab some
    dark frames. Then move lmir_fw2 back to where it was.'''
    pi.setINDI("LBTO.OffsetPointing.CoordSys", coords,
           "LBTO.OffsetPointing.OffsetX", DecOff, 
           "LBTO.OffsetPointing.OffsetY", RAOff, 
           "LBTO.OffsetPointing.Side", side, 
           "LBTO.OffsetPointing.Type", absrel, 
           timeout=400,
           wait=False) 
	###where's FW2
    lightIntParams = pi.getINDI("LMIRCAM.acquire")
    currentFlag = pi.getINDI("LMIRCAM.EditFITS")
    setFLAG('DRK')
    fw2_pos = pi.getINDI("Lmir.lmir_FW2_status.PosNum", wait=True)
    pi.setINDI("Lmir.lmir_FW2.command", int(fw2_pos)+25000, timeout=20, wait=True)
    get_lmircam_frames(*darkIntParams,use_bg=True)
    pi.setINDI("Lmir.lmir_FW2.command", int(fw2_pos), timeout=20, wait=True)
    pi.setINDI(currentFlag)
    wait4AORunning(side)

def get_ales_frames(dit, coadds, nseqs, side='left', use_bg=0, save_data=True, flag=''):           #  Standard function to execute LMIRCam exosures.
    pi.setINDI("LMIRCAM.stop.value=Off")                                        #  Make sure the camera is not running.
    
    if save_data:                                                               #  If data saving is enabled.
        pi.setINDI("LMIRCAM.enable_save.value=On")                              #  Turn on data saving.
    else:                                                                       #  Else.
        pi.setINDI("LMIRCAM.enable_save.value=Off")                             #  Turn off data saving.
    
    bxcur, bycur, bxreq, byreq = get_baysidestage_positions(side)
    print 'bayside stage positions, xcur, ycur, xreq, yreq:'
    print bxcur, bycur, bxreq, byreq
    pi.setINDI("LMIRCAM.EditFITS.Keyword=FLAG;Value=%s ;Comment=observation type"%flag, wait=False)
    pi.setINDI("LMIRCAM.EditFITS.Keyword=BAYXCUR;Value=%f ;Comment=current baysidex position"%bxcur,wait=True)
    pi.setINDI("LMIRCAM.EditFITS.Keyword=BAYXREQ;Value=%f ;Comment=requested baysidex position"%bxreq,wait=True)
    pi.setINDI("LMIRCAM.EditFITS.Keyword=BAYYCUR;Value=%f ;Comment=current baysidex position"%bycur,wait=True)
    pi.setINDI("LMIRCAM.EditFITS.Keyword=BAYYREQ;Value=%f ;Comment=requested baysidex position"%byreq,wait=True)
    #  Now run the camera with the parameters provided: print dit, coadds, nseqs, use_bg
    pi.setINDI("LMIRCAM.acquire.int_time=%f;num_coadds=%i;num_seqs=%i;enable_bg=%i;is_bg=0;is_cont=0"%(dit, coadds, nseqs, use_bg), timeout = 2.0 * dit * nseqs)

def get_baysidestage_positions(side):
    sidestrs = {'left':('sx','L'),'right':('dx','R')}[side]
    bayxcur = pi.getINDI("indi_"+sidestrs[0]+"_wfs_msgd.baysidex___"+sidestrs[1]+"___POS___CUR.value")
    bayxreq = pi.getINDI("indi_"+sidestrs[0]+"_wfs_msgd.baysidex___"+sidestrs[1]+"___POS___REQ.value")
    bayycur = pi.getINDI("indi_"+sidestrs[0]+"_wfs_msgd.baysidey___"+sidestrs[1]+"___POS___CUR.value")
    bayyreq = pi.getINDI("indi_"+sidestrs[0]+"_wfs_msgd.baysidey___"+sidestrs[1]+"___POS___REQ.value")
    return bayxcur, bayycur, bayxreq, bayyreq
