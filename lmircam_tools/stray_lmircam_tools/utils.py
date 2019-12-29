from lmircam_tools import *

def rawbg():
    pi.setINDI("LMIRCAM.Command.text","rawbg")

def setFLAG(string):
    pi.setINDI("LMIRCAM.EditFITS.Keyword=FLAG;Value=%s;Comment=Script Note"%string, wait=False)

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

def nod(XOff,YOff,side,Coords="DETXY",Type="REL"):
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
