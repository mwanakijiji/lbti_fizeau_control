# This reads in csvs of injected and retrieved OPD, tip, and tilt in Fizeau PSFs

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

wavel_um = 3.87
PS = 10.7

def plot_analysis(csv_file):
    '''
    INPUTS:
    csv_file: the csv containing the results of the analysis
    '''

    # read in all the data (without choosing which to plot quite yet)
    # OPD variations alone
    df = pd.read_csv(csv_file,
                    names=["file_name","opd_inject",
                        "tip_inject","tilt_inject",
                        "x_shift_inject","y_shift_inject",
                        "opd_retrieve","tip_retrieve",
                        "tilt_retrieve","x_shift_retrieve",
                        "y_shift_retrieve"])
    df_sort = df.sort_values(by="file_name").reset_index()

    # get quantities onto common convention

    # N.b.: the retrieved OPD is with respect to the nearest wavelength-multiple
    # of OPD

    # injected quantities
    opd_inj_unwrapped = df_sort["opd_inject"]
    opd_inj_simple_mod = np.mod(df_sort["opd_inject"],wavel_um)
    opd_inj_positive_side = np.copy(opd_inj_simple_mod)
    opd_inj_negative_side = np.copy(opd_inj_simple_mod)
    tip_inj_unwrapped = df_sort["tip_inject"]
    tip_inj_simple_mod = np.mod(df_sort["tip_inject"],2*PS)
    tip_inj_positive_side = np.copy(tip_inj_simple_mod)
    tip_inj_negative_side = np.copy(tip_inj_simple_mod)
    tilt_inj_unwrapped = df_sort["tilt_inject"]
    tilt_inj_simple_mod = np.mod(df_sort["tilt_inject"],2*PS)
    tilt_inj_positive_side = np.copy(tilt_inj_simple_mod)
    tilt_inj_negative_side = np.copy(tilt_inj_simple_mod)
    delx_inj = df_sort["x_shift_inject"]
    dely_inj = df_sort["y_shift_inject"]

    # retrieved quantities
    opd_ret = df_sort["opd_retrieve"]
    tip_ret = df_sort["tip_retrieve"]
    tilt_ret = df_sort["tilt_retrieve"]
    delx_ret = df_sort["x_shift_retrieve"]
    dely_ret = df_sort["y_shift_retrieve"]

    # injected OPDs, after wrapping, and within an OPD range of (0, +lambda/2)
    opd_inj_positive_side[
        ~np.logical_and(opd_inj_simple_mod >= 0,opd_inj_simple_mod < wavel_um/2)] = np.nan
    # injected OPDs, after wrapping, and within an OPD range of (-lambda/2, 0)
    opd_inj_negative_side[
        ~np.logical_and(opd_inj_simple_mod >= wavel_um/2,opd_inj_simple_mod < wavel_um)] = np.nan
    # translate the 'negative' side OPDs
    opd_inj_negative_side -= wavel_um

    # injected tip, after wrapping, and within a tilt range of (0, +PS)
    tip_inj_positive_side[
        ~np.logical_and(tip_inj_simple_mod >= 0,tip_inj_simple_mod < PS)] = np.nan
    # injected tips, after wrapping, and within a tip range of (-PS, 0)
    tip_inj_negative_side[
        ~np.logical_and(tip_inj_simple_mod >= PS,tip_inj_simple_mod < 2*PS)] = np.nan
    # translate the 'negative' side tips
    tip_inj_negative_side -= 2*PS

    # injected tilt, after wrapping, and within a tilt range of (0, +PS)
    tilt_inj_positive_side[
        ~np.logical_and(tilt_inj_simple_mod >= 0,tilt_inj_simple_mod < PS)] = np.nan
    # injected tilts, after wrapping, and within a tip range of (-PS, 0)
    tilt_inj_negative_side[
        ~np.logical_and(tilt_inj_simple_mod >= PS,tilt_inj_simple_mod < 2*PS)] = np.nan
    # translate the 'negative' side tilts
    tilt_inj_negative_side -= 2*PS

    #############################################
    ###### BEGIN PLOT OF WRAPPED QUANTITIES #####
    #import ipdb; ipdb.set_trace()
    fig = plt.figure(figsize=(12,10), constrained_layout=True)
    gs = gridspec.GridSpec(6, 1, figure=fig)

    #import ipdb; ipdb.set_trace()
    # OPD
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(opd_inj_positive_side, color = "b", linewidth=10, alpha=0.25, label="Injected")
    ax1.plot(opd_inj_negative_side, color = "b", linewidth=10, alpha=0.25)
    ax1.plot(opd_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1, label="Retrieved")
    ax1.axhline(y=wavel_um/2, linestyle = "--")
    ax1.axhline(y=-wavel_um/2, linestyle = "--")
    ax1.annotate("$+\lambda/2$", xy=(0, 0), xytext=(40, 2.1), textcoords="data", fontsize = 15.0)
    ax1.annotate("$-\lambda/2$", xy=(0, 0), xytext=(40, -2.3), textcoords="data", fontsize = 15.0)
    #ax1.set_title("OPD, wrapped\n(blue = injected; red = retrieved)", fontsize = 18.0)
    #ax1.set_xlim([0,10000])
    #ax1.set_ylim([-2.5,2.5])
    #ax1.set_ylim([np.nanmin(opd_inj_unwrapped),np.nanmax(opd_inj_unwrapped)])
    ax1.set_ylabel("OPD ($\mu$m)", fontsize = 15.0)
    ax1.set_xlabel("Elapsed time (sec)", fontsize = 15.0)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.tick_params(labelsize=14)
    ax1.legend()

    # tip
    ax2 = fig.add_subplot(gs[1, :])
    #ax2.plot(elapsed_time_inj, tip_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax2.plot(tip_inj_positive_side, color = "b", linewidth=10, alpha=0.25)
    ax2.plot(tip_inj_negative_side, color = "b", linewidth=10, alpha=0.25)
    ax2.plot(tip_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax2.axhline(y=10.7, linestyle = "--")
    ax2.axhline(y=-10.7, linestyle = "--")
    ax2.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax2.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax2.set_title("Tip, wrapped (y)", fontsize = 18.0)
    #ax2.set_xlim([0,10000])
    #ax2.set_ylim([-14,14])
    #ax2.set_ylim([2*np.nanmin(tip_inj_unwrapped),2*np.nanmax(tip_inj_unwrapped)])
    ax2.set_ylabel("Tip (y),\nwrapped  (mas)", fontsize = 15.0)
    ax2.set_xlabel("Elapsed time (sec)", fontsize = 15.0)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.tick_params(labelsize=14)

    # tilt
    ax3 = fig.add_subplot(gs[2, :])
    #ax3.plot(elapsed_time_inj, tilt_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax3.plot(tilt_inj_positive_side, color = "b", linewidth=10, alpha=0.25)
    ax3.plot(tilt_inj_negative_side, color = "b", linewidth=10, alpha=0.25)
    ax3.plot(tilt_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax3.axhline(y=10.7, linestyle = "--")
    ax3.axhline(y=-10.7, linestyle = "--")
    ax3.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax3.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax3.set_title("Tilt, wrapped (x)", fontsize = 18.0)
    #ax3.set_xlim([0,10000])
    #ax3.set_ylim([-14,14])
    #ax3.set_ylim([2*np.nanmin(tilt_inj_unwrapped),2*np.nanmax(tilt_inj_unwrapped)])
    ax3.set_ylabel("Tilt (x),\nwrapped (mas)", fontsize = 15.0)
    ax3.axes.get_xaxis().set_visible(False)
    ax3.tick_params(labelsize=14)

    # x
    ax4 = fig.add_subplot(gs[3, :])
    #ax4.plot(elapsed_time_inj, tilt_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax4.plot(delx_inj, color = "b", linewidth=10, alpha=0.25)
    ax4.plot(delx_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax4.axhline(y=10.7, linestyle = "--")
    ax4.axhline(y=-10.7, linestyle = "--")
    ax4.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax4.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax4.set_title("Tilt, wrapped (x)", fontsize = 18.0)
    #ax4.set_xlim([0,10000])
    #ax4.set_ylim([-14,14])
    #ax4.set_ylim([2*np.nanmin(delx_inj),2*np.nanmax(delx_inj)])
    ax4.set_ylabel("Translation ($\Delta$x),\n(pix)", fontsize = 15.0)
    ax4.axes.get_xaxis().set_visible(False)
    ax4.tick_params(labelsize=14)

    # y
    ax5 = fig.add_subplot(gs[4, :])
    #ax5.plot(elapsed_time_inj, tilt_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax5.plot(dely_inj, color = "b", linewidth=10, alpha=0.25)
    ax5.plot(dely_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax5.axhline(y=10.7, linestyle = "--")
    ax5.axhline(y=-10.7, linestyle = "--")
    ax5.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax5.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax5.set_title("Tilt, wrapped (x)", fontsize = 18.0)
    #ax5.set_xlim([0,10000])
    #ax5.set_ylim([-14,14])
    #ax5.set_ylim([2*np.nanmin(dely_inj),2*np.nanmax(dely_inj)])
    ax5.set_ylabel("Translation ($\Delta$y),\n(pix)", fontsize = 15.0)
    ax5.set_xlabel("Consecutive frame num (unitless)", fontsize = 15.0)
    ax5.tick_params(labelsize=14)

    plt.suptitle(csv_file)

    plt.savefig("test_wrapped.pdf")
    plt.clf()
    plt.close()
    ##### END PLOT OF WRAPPED QUANTITIES #####
    ##########################################################

    ##########################################################
    ###### BEGIN PLOT OF ABSOLUTE / UNWRAPPED QUANTITIES #####

    fig = plt.figure(figsize=(12,10), constrained_layout=True)
    gs = gridspec.GridSpec(6, 1, figure=fig)

    # OPD
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(opd_inj_unwrapped, color = "b", linewidth=10, alpha=0.25, label="Injected")
    ax1.plot(opd_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1, label="Retrieved")
    ax1.axhline(y=wavel_um/2, linestyle = "--")
    ax1.axhline(y=-wavel_um/2, linestyle = "--")
    ax1.annotate("$+\lambda/2$", xy=(0, 0), xytext=(40, 2.1), textcoords="data", fontsize = 15.0)
    ax1.annotate("$-\lambda/2$", xy=(0, 0), xytext=(40, -2.3), textcoords="data", fontsize = 15.0)
    #ax1.set_title("OPD, wrapped\n(blue = injected; red = retrieved)", fontsize = 18.0)
    #ax1.set_xlim([0,10000])
    #ax1.set_ylim([2*np.nanmin(opd_inj_unwrapped),2*np.nanmax(opd_inj_unwrapped)])
    ax1.set_ylabel("OPD ($\mu$m)", fontsize = 15.0)
    ax1.set_xlabel("Elapsed time (sec)", fontsize = 15.0)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.tick_params(labelsize=14)
    ax1.legend()

    # tip
    ax2 = fig.add_subplot(gs[1, :])
    #ax2.plot(elapsed_time_inj, tip_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax2.plot(tip_inj_unwrapped, color = "b", linewidth=10, alpha=0.25)
    ax2.plot(tip_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax2.axhline(y=10.7, linestyle = "--")
    ax2.axhline(y=-10.7, linestyle = "--")
    ax2.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax2.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax2.set_title("Tip, wrapped (y)", fontsize = 18.0)
    #ax2.set_xlim([0,10000])
    #ax2.set_ylim([2*np.nanmin(tip_inj_unwrapped),2*np.nanmax(tip_inj_unwrapped)])
    ax2.set_ylabel("Tip (y),\nabsolute  (mas)", fontsize = 15.0)
    ax2.set_xlabel("Elapsed time (sec)", fontsize = 15.0)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.tick_params(labelsize=14)

    # tilt
    ax3 = fig.add_subplot(gs[2, :])
    #ax3.plot(elapsed_time_inj, tilt_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax3.plot(tilt_inj_unwrapped, color = "b", linewidth=10, alpha=0.25)
    ax3.plot(tilt_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax3.axhline(y=10.7, linestyle = "--")
    ax3.axhline(y=-10.7, linestyle = "--")
    ax3.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax3.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax3.set_title("Tilt, wrapped (x)", fontsize = 18.0)
    #ax3.set_xlim([0,10000])
    #ax3.set_ylim([2*np.nanmin(tilt_inj_unwrapped),2*np.nanmax(tilt_inj_unwrapped)])
    ax3.set_ylabel("Tilt (x),\nabsolute (mas)", fontsize = 15.0)
    ax3.axes.get_xaxis().set_visible(False)
    ax3.tick_params(labelsize=14)

    # x
    ax4 = fig.add_subplot(gs[3, :])
    #ax4.plot(elapsed_time_inj, tilt_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax4.plot(delx_inj, color = "b", linewidth=10, alpha=0.25)
    ax4.plot(delx_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax4.axhline(y=10.7, linestyle = "--")
    ax4.axhline(y=-10.7, linestyle = "--")
    ax4.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax4.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax4.set_title("Tilt, wrapped (x)", fontsize = 18.0)
    #ax4.set_xlim([0,10000])
    #ax4.set_ylim([2*np.nanmin(delx_inj),2*np.nanmax(delx_inj)])
    ax4.set_ylabel("Translation ($\Delta$x),\n(pix)", fontsize = 15.0)
    ax4.axes.get_xaxis().set_visible(False)
    ax4.tick_params(labelsize=14)

    # y
    ax5 = fig.add_subplot(gs[4, :])
    #ax5.plot(elapsed_time_inj, tilt_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax5.plot(dely_inj, color = "b", linewidth=10, alpha=0.25)
    ax5.plot(dely_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax5.axhline(y=10.7, linestyle = "--")
    ax5.axhline(y=-10.7, linestyle = "--")
    ax5.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax5.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax5.set_title("Tilt, wrapped (x)", fontsize = 18.0)
    #ax5.set_xlim([0,10000])
    #ax5.set_ylim([2*np.nanmin(dely_inj),2*np.nanmax(dely_inj)])
    ax5.set_ylabel("Translation ($\Delta$y),\n(pix)", fontsize = 15.0)
    ax5.set_xlabel("Consecutive frame num (unitless)", fontsize = 15.0)
    ax5.tick_params(labelsize=14)

    plt.suptitle(csv_file)

    plt.savefig("test_absolute.pdf")
