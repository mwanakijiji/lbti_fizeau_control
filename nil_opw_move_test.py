from pyindi import *
pi = PyINDI(verbose=False)

movement_dir = "pupils_in"

if (movement_dir == "pupils_in"):  # I dont know if this is left or right on LMIR!
    print("Moving the NIL OPW pupils in")
    pi.setINDI("Warm.NIL_OPW_status.PosName", "NomPupilsIn", wait=True)
elif (movement_dir == "pupils_out"):
    print("Moving the NIL OPW pupils out")
    pi.setINDI("Warm.NIL_OPW_status.PosName", "NomPupilsOut", wait=True)
else:
    print("No legit choice of NIL OPW pupil!!")

print("Now move the telescopes either left or right 5 asec; edit script to get direction right")
