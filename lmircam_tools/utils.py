from lmircam_tools import *
import os

#def rawbg():
#    pi.setINDI("LMIRCAM.Command.text","rawbg")
 
def setFLAG(string):
    pi.setINDI("lmircam_save.EditFITS.Keyword=FLAG;Value=%s;Comment=Script Note"%string, wait=False)

def wait4AORunning(side):
    pistr={ 'left':["LBTO.Dictionary.Name=L_AOStatus;Value="],
            'right':["LBTO.Dictionary.Name=R_AOStatus;Value="],
            'both':["LBTO.Dictionary.Name=L_AOStatus;Value=",
                    "LBTO.Dictionary.Name=R_AOStatus;Value="]}[side]
    while True:
        status=[]
        for side_string in pistr:
            pi.setINDI(side_string)
            state = pi.getINDI("LBTO.Dictionary.Value")
            status.append(state)
        time.sleep(0.02)
        if all([s == "AORunning" for s in status]):
            break

####  The following function is a workaround version of the 'nod' function to fix unstable
####  DX nods.  Remove/rename once the issue is fixed and rename the 'original' nod function
####  below (now called 'nod_regular') to 'nod'.

def nod(XOff,YOff,side,Coords="DETXY",Type="REL"):
    '''side should be one of 'left', 'right', or 'both' '''
####  Here are some lines to fix DX issue.
    if side == 'right' or side == 'both':    # This is only relevant if DX is involved.
      low_gain_file = 'gain0.1_10modi.fits'  # File name for low gain file.
      settling_time = 3                      # Time the AO should be closed in low gains before setting high gain again.
      full_gain_path = pi.getINDI('indi_dx_adsec_msgd.ADSEC___R___G_GAIN_A.value')  # Path of the full gain file.
      low_gain_path = os.path.join(os.path.dirname(full_gain_path), low_gain_file)  # Path of the low gain file (same directory, but above file name).
      print 'Setting 10-mode gain.'          # Some status info.
      pi.setINDI('ao_indi_dx_ice.load_gain_tag_a.tag', low_gain_path, wait=True)    # Set low gains before nodding.
####  Those were some lines to fix DX issue.
    print 'nodding %f x, %f y on side %s' % (XOff, YOff, side)
    pi.setINDI("LBTO.OffsetPointing.CoordSys", Coords,
               "LBTO.OffsetPointing.OffsetX", XOff,
               "LBTO.OffsetPointing.OffsetY", YOff,
               "LBTO.OffsetPointing.Side", side,
               "LBTO.OffsetPointing.Type", Type,
               timeout=400,
               wait=True)
    wait4AORunning(side)
####  Here are some more lines to fix DX issue.
    if side == 'right' or side == 'both':  # This is only relevant if DX is involved.
      print 'Waiting a few seconds...'     # Some status info.
      sleep(settling_time)                 # Sleep to allow AO to settle in low gain.
      pi.setINDI('ao_indi_dx_ice.load_gain_tag_a.tag', full_gain_path, wait=True)  # Set gains back to original, high values.
      print 'Ok to continue.'              # All done.
####  Those were some more lines to fix DX issue.


####  The following is the 'normal' nod without the workaround to fix unstable DX nods.
####  Roll back once the issue is fixed by renaming the function 'nod' above to something
####  else and this one to 'nod'.

def nod_regular(XOff,YOff,side,Coords="DETXY",Type="REL"):
    '''side should be one of 'left', 'right', or 'both' '''
    print 'nodding %f x, %f y on side %s' % (XOff, YOff, side)
    pi.setINDI("LBTO.OffsetPointing.CoordSys", Coords,
               "LBTO.OffsetPointing.OffsetX", XOff,
               "LBTO.OffsetPointing.OffsetY", YOff,
               "LBTO.OffsetPointing.Side", side,
               "LBTO.OffsetPointing.Type", Type,
               timeout=400,
               wait=True)
    wait4AORunning(side)

