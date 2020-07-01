# This reads in a bunch of MTFs and PTFs and plots stuff to try to diagnose
# degeneracies between different aberrations
# (note this assumes 500x500 cutouts for now; 2020 June 29, E.S.)

from modules import *
from modules import plot_injected_retrieved
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import multiprocessing as mp
import numpy as np
import pandas as pd
import pyfits
import glob
import os
import time

# subsets for testing
#path_stem = "/vol_c/synthetic_fizeau/otfs_trial1_opd_tip_tilt/"
#path_stem = "/vol_c/synthetic_fizeau/otfs_trial2_opd/"
#path_stem = "/vol_c/synthetic_fizeau/otfs_trial3_tip/"
#path_stem = "/vol_c/synthetic_fizeau/otfs_trial4_tilt/"
#path_stem = "/vol_c/synthetic_fizeau/otfs_trial5_yx/"
#path_stem = "/vol_c/synthetic_fizeau/otfs_trial6_opd_tip_tilt_yx/"
path_stem = "/vol_c/synthetic_fizeau/otfs_all_subsets/"
path_for_plots = path_stem + "plots/"


def worker(file_name, q):
    # reads in a single FITS file, calculates some quantities, and returns them

    string_basename = os.path.basename(file_name)
    plot_write_name = path_for_plots + "png_xsecs_" + string_basename.split(".")[0] + ".pdf"

    #print("Making plot for " + string_basename)

    f = pyfits.open(file_name)
    img = f[0].data
    header = f[0].header # contains injected values

    '''
    opd_inj_unwrapped = header["OPD_UM"]
    inj_tip = header["TIPY_MAS"]
    inj_tilt = header["TILTXMAS"]
    if np.logical_and(header["X_SHIFT"],header["Y_SHIFT"]):
        inj_delx = header["X_SHIFT"]
        inj_dely = header["Y_SHIFT"]
    else:
        inj_delx = 0.
        inj_dely = 0.
    '''

    # take x-sec in x
    x_sec = img[250,:]

    # take x-sec in y across each lobe
    y_sec_left = img[:,153]
    y_sec_middle = img[:,250]
    y_sec_right = img[:,348]

    # calculate wrapped OPD
    # injected quantities
    '''
    opd_inj_simple_mod = np.mod(opd_inj_unwrapped,wavel_um)
    opd_inj_positive_side = np.copy(opd_inj_simple_mod)
    opd_inj_negative_side = np.copy(opd_inj_simple_mod)

    # injected OPDs, after wrapping, and within an OPD range of (0, +lambda/2)
    opd_inj_positive_side[
        ~np.logical_and(opd_inj_simple_mod >= 0,opd_inj_simple_mod < wavel_um/2)] = np.nan
    # injected OPDs, after wrapping, and within an OPD range of (-lambda/2, 0)
    opd_inj_negative_side[
        ~np.logical_and(opd_inj_simple_mod >= wavel_um/2,opd_inj_simple_mod < wavel_um)] = np.nan
    # translate the 'negative' side OPDs
    opd_inj_negative_side -= wavel_um
    '''

    # plot it: x-sec in x across the top, and 3 x-sec in y along the bottom
    fig=plt.figure(figsize=(9,8))
    gs=GridSpec(2,3) # 2 rows, 3 columns
    ax1=fig.add_subplot(gs[0,:]) # Second row, span all columns
    ax2=fig.add_subplot(gs[1,0]) # First row, first column
    ax3=fig.add_subplot(gs[1,1]) # First row, second column
    ax4=fig.add_subplot(gs[1,2]) # First row, third column

    ax1.plot(x_sec)
    ax2.plot(y_sec_left)
    ax3.plot(y_sec_middle)
    ax4.plot(y_sec_right)

    ax1.set_title("x-sec in x")
    ax1.set_xlabel("pix$_{DFT}$ (y)")
    ax1.set_ylabel("MTF or PTF")
    ax2.set_title("x-sec in y, left")
    ax2.set_xlabel("pix$_{DFT}$ (x)")
    ax2.set_ylabel("MTF or PTF")
    ax3.set_title("x-sec in y, middle")
    ax3.set_xlabel("pix$_{DFT}$ (x)")
    ax4.set_title("x-sec in y, right")
    ax4.set_xlabel("pix$_{DFT}$ (x)")
    plt.suptitle(os.path.basename(file_name))
    '''
    plt.suptitle(os.path.basename(file_name) + \
        "\ninjected tip: " + str(inj_tip) + \
        "\ninjected tilt: " + str(inj_tilt) + \
        "\ninjected opd: " + str(opd_inj_unwrapped) + \
        "\ninjected del_x: " + str(inj_delx) + \
        "\ninjected del_y: " + str(inj_dely))
    '''

    plt.tight_layout()

    fig.savefig(plot_write_name)

    return

def listener(q):
    '''listens for messages on the q, writes to file. '''

    with open(fn_write, 'w') as f:
        while 1:
            m = q.get()
            if m == 'kill':
                break
            f.write(str(m) + '\n')
            f.flush()

def main():
    '''
    Grab frames from a directory and retrieve the aberrations in the PSF
    '''

    # get list of file names together
    files_list = glob.glob(path_stem + "*.fits")
    #print("Files to analyze: ")
    #print(files_list)

    #must use Manager queue here, or will not work
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() + 2)

    #put listener to work first
    watcher = pool.apply_async(listener, (q,))

    #fire off workers
    jobs = []
    for i in range(len(files_list)):
        file_name_fits = files_list[i]
        job = pool.apply_async(worker, (file_name_fits, q))
        jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    #now we are done, kill the listener
    q.put('kill')
    pool.close()
    pool.join()


if __name__ == "__main__":
    main() # analysis of frames
