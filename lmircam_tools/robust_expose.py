from lmircam_tools import *
from multiprocessing import Process
from print_tools import *

class LMIRCamException(Exception):
    pass

def snoop_for_badrows(inner_loop_max, timeout=900):
    n_bads = 0
    print 'snoop: if more than %i are bad will raise an exception' % int(0.2*inner_loop_max)
    for ii in xrange(inner_loop_max):
        f = pi.getFITS("LMIRCAM.DSImage.File", timeout=0)
        print ii, f[0].header['FILENAME']
        if 'BADROWS' in f[0].header.keys():
            print 'found one'
            n_bads +=1
            last_few +=1
            print "n found: ", n_bads
        else: 
            last_few = 0
        if n_bads > int(0.2*inner_loop_max):
            pi.setINDI("LMIRCAM.CommandOOB.text","!ALLSTOP")
            raise LMIRCamException('Too many bad rows total')
        if last_few > 5:
            pi.setINDI("LMIRCAM.CommandOOB.text","!ALLSTOP")
            raise LMIRCamException('Too many bad rows in sequence')

def get_lmircam_frames_robust(dit, coadds, nseqs, 
                              inner_loop_max=10, 
                              save_data=True):
    dit = dit
    coadds = int(coadds)
    nseqs = int(nseqs)
    inner_loop_max = int(inner_loop_max)
    N_outer_loops = nseqs / inner_loop_max
    remainder = nseqs % inner_loop_max
    lbtintparams = (dit, coadds, inner_loop_max)

    pi.setINDI("LMIRCAM.Command.text", "0 contacq")
    if save_data:
        pi.setINDI("LMIRCAM.Command.text", "1 savedata")
    else:
        pi.setINDI("LMIRCAM.Command.text", "0 savedata")
        info('Data saving is turned OFF')
    pi.setINDI("LMIRCAM.Command.text", "%f %i %i lbtintpar" % lbtintparams, wait=True)

    for ii in xrange(N_outer_loops):
        try: 
            snoop_bads = Process(target=snoop_for_badrows,args=(inner_loop_max,))
            snoop_bads.start()
            pi.setINDI("LMIRCAM.Command.text","go",timeout=(8*dit*coadds*inner_loop_max),wait=True)
            print 'Exit code before: ', snoop_bads.exitcode
            snoop_bads.terminate()
            print 'Exit code after: ', snoop_bads.exitcode
        except LMIRCamException:
            info("LMIRCam Exception caught")
            raw_input()
            
    if remainder > 0:
        pi.setINDI("LMIRCAM.Command.text", "%f %i %i lbtintpar" % (dit, coadds, remainder), wait=True)
        snoop_bads = Process(target=snoop_for_badrows,args=(inner_loop_max,))
        snoop_bads.start()
        pi.setINDI("LMIRCAM.Command.text","go",timeout=(8*dit*coadds*nseqs),wait=True) 
        print 'Exit code before: ', snoop_bads.exitcode
        snoop_bads.terminate()
        print 'Exit code after: ', snoop_bads.exitcode


