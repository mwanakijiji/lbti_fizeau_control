from pyindi import *
from lmircam_tools.utils import nod

pi = PyINDI(verbose=False)

movement_dir = "pupils_right"

if (movement_dir == "pupils_right"):  # I dont know if this is left or right on LMIR!
    print("Moving the NIL OPW pupils in")
    #pi.setINDI("Warm.NIL_OPW.command", "NomPupilsIn", wait=True)
    pi.setINDI("Warm.NIL_OPW.command", "-416700", wait=True) # goes to right on Phasecam
elif (movement_dir == "pupils_left"):
    print("Moving the NIL OPW pupils out")
    #pi.setINDI("Warm.NIL_OPW.command", "NomPupilsOut", wait=True)
    pi.setINDI("Warm.NIL_OPW.command", "0", wait=True) # goes to left on Phasecam
else:
    print("No legit choice of NIL OPW pupil!!")

print("Now move the telescopes either left or right 5 asec; edit script to get direction right")

## HERE ARE THE TELESCOPE MOVE COMMANDS
#nod_x = 5 (or -5)
#nod_y = 0
#side = "both"
#nod(nod_x, nod_y, side)
