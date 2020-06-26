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

    # retrieved quantities
    opd_ret = df_sort["opd_retrieve"]
    tip_ret = df_sort["tip_retrieve"]
    tilt_ret = df_sort["tilt_retrieve"]

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


    fig = plt.figure(figsize=(12,10), constrained_layout=True)
    gs = gridspec.GridSpec(6, 1, figure=fig)


    ### CHOOSE FIRST PANEL FOR THE UNWRAPPED QUANTITY OF INTEREST
    # injected OPD, unwrapped
    '''
    ax0 = fig.add_subplot(gs[0, :])
    ax0.plot(elapsed_time_inj, opd_inj_unwrapped, color = "b", linewidth=10, alpha=0.15)
    ax0.axhline(y=wavel_um/2, linestyle = "--")
    ax0.axhline(y=-wavel_um/2, linestyle = "--")
    ax0.annotate("$+\lambda/2$", xy=(0, 0), xytext=(40, 2.1), textcoords="data")
    ax0.annotate("$-\lambda/2$", xy=(0, 0), xytext=(40, -2.3), textcoords="data")
    ax0.set_title("Unwrapped injected OPD")
    ax0.set_xlim([0,10000])
    #ax0.set_ylim([-2.5,2.5])
    ax0.set_ylabel("OPD ($\mu$m)")
    ax0.set_xlabel("Elapsed time (sec)")
    '''

    # injected tip, unwrapped
    '''
    ax0 = fig.add_subplot(gs[0, :])
    ax0.plot(elapsed_time_inj, tip_inj_unwrapped, color = "b", linewidth=10, alpha=0.15)
    ax0.axhline(y=10.7, linestyle = "--")
    ax0.axhline(y=-10.7, linestyle = "--")
    ax0.annotate("$+PS$", xy=(0, 0), xytext=(40, 12), textcoords="data")
    ax0.annotate("$-PS$", xy=(0, 0), xytext=(40, -18), textcoords="data")
    ax0.set_title("Unwrapped injected tip")
    ax0.set_xlim([0,10000])
    ax0.set_ylim([-20,60])
    ax0.set_ylabel("Tip (mas)")
    ax0.set_xlabel("Elapsed time (sec)")
    '''

    # injected tilt, unwrapped
    ax0 = fig.add_subplot(gs[0, :])
    #ax0.axvline(x=400, linestyle=":", color="gray")
    #ax0.axvline(x=1200, linestyle=":", color="gray")
    ax0.plot(tilt_inj_unwrapped, color = "b", linewidth=10, alpha=0.15, label="Injected")
    ax0.plot(opd_ret-1000, color="red", marker='o', linestyle='None', markersize=2, alpha=1, label="Retrieved") # fake; just for legend
    ax0.axhline(y=10.7, linestyle = "--")
    ax0.axhline(y=-10.7, linestyle = "--")
    ax0.annotate("$+PS$", xy=(0, 0), xytext=(40, 12), textcoords="data", fontsize = 15.0)
    ax0.annotate("$-PS$", xy=(0, 0), xytext=(40, -18), textcoords="data", fontsize = 15.0)
    #ax0.set_title("Unwrapped injected tilt", fontsize = 18.0)
    ax0.set_xlim([0,10000])
    ax0.set_ylim([-20,40])
    ax0.set_ylabel("Injected tilt,\nabsolute (mas)", fontsize = 15.0)
    ax0.set_xlabel("Elapsed time (sec)", fontsize = 15.0)
    ax0.axes.get_xaxis().set_visible(False)
    ax0.tick_params(labelsize=14)
    ax0.legend()


    ### END FIRST PANEL

    # OPD
    ax1 = fig.add_subplot(gs[1, :])
    ax1.axvline(x=400, linestyle=":", color="gray")
    ax1.axvline(x=1200, linestyle=":", color="gray")
    ax1.plot(opd_inj_positive_side, color = "b", linewidth=10, alpha=0.25)
    ax1.plot(opd_inj_negative_side, color = "b", linewidth=10, alpha=0.25)
    ax1.plot(opd_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax1.axhline(y=wavel_um/2, linestyle = "--")
    ax1.axhline(y=-wavel_um/2, linestyle = "--")
    ax1.annotate("$+\lambda/2$", xy=(0, 0), xytext=(40, 2.1), textcoords="data", fontsize = 15.0)
    ax1.annotate("$-\lambda/2$", xy=(0, 0), xytext=(40, -2.3), textcoords="data", fontsize = 15.0)
    #ax1.set_title("OPD, wrapped\n(blue = injected; red = retrieved)", fontsize = 18.0)
    ax1.set_xlim([0,10000])
    ax1.set_ylim([-2.5,2.5])
    ax1.set_ylabel("OPD ($\mu$m)", fontsize = 15.0)
    ax1.set_xlabel("Elapsed time (sec)", fontsize = 15.0)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.tick_params(labelsize=14)

    # tip
    ax2 = fig.add_subplot(gs[2, :])
    ax2.axvline(x=400, linestyle=":", color="gray")
    ax2.axvline(x=1200, linestyle=":", color="gray")
    #ax2.plot(elapsed_time_inj, tip_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax2.plot(tip_inj_positive_side, color = "b", linewidth=10, alpha=0.25)
    ax2.plot(tip_inj_negative_side, color = "b", linewidth=10, alpha=0.25)
    ax2.plot(tip_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax2.axhline(y=10.7, linestyle = "--")
    ax2.axhline(y=-10.7, linestyle = "--")
    ax2.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax2.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax2.set_title("Tip, wrapped (y)", fontsize = 18.0)
    ax2.set_xlim([0,10000])
    ax2.set_ylim([-14,14])
    ax2.set_ylabel("Tip (y),\nwrapped  (mas)", fontsize = 15.0)
    ax2.set_xlabel("Elapsed time (sec)", fontsize = 15.0)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.tick_params(labelsize=14)

    # tilt
    ax3 = fig.add_subplot(gs[3, :])
    ax3.axvline(x=400, linestyle=":", color="gray")
    ax3.axvline(x=1200, linestyle=":", color="gray")
    #ax3.plot(elapsed_time_inj, tilt_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax3.plot(tilt_inj_positive_side, color = "b", linewidth=10, alpha=0.25)
    ax3.plot(tilt_inj_negative_side, color = "b", linewidth=10, alpha=0.25)
    ax3.plot(tilt_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax3.axhline(y=10.7, linestyle = "--")
    ax3.axhline(y=-10.7, linestyle = "--")
    ax3.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax3.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax3.set_title("Tilt, wrapped (x)", fontsize = 18.0)
    ax3.set_xlim([0,10000])
    ax3.set_ylim([-14,14])
    ax3.set_ylabel("Tilt (x),\nwrapped (mas)", fontsize = 15.0)
    ax3.axes.get_xaxis().set_visible(False)
    ax3.tick_params(labelsize=14)

    # x
    ax4 = fig.add_subplot(gs[4, :])
    ax4.axvline(x=400, linestyle=":", color="gray")
    ax4.axvline(x=1200, linestyle=":", color="gray")
    #ax4.plot(elapsed_time_inj, tilt_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax4.plot(tilt_inj_positive_side, color = "b", linewidth=10, alpha=0.25)
    ax4.plot(tilt_inj_negative_side, color = "b", linewidth=10, alpha=0.25)
    ax4.plot(tilt_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax4.axhline(y=10.7, linestyle = "--")
    ax4.axhline(y=-10.7, linestyle = "--")
    ax4.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax4.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax4.set_title("Tilt, wrapped (x)", fontsize = 18.0)
    ax4.set_xlim([0,10000])
    ax4.set_ylim([-14,14])
    ax4.set_ylabel("Translation ($\Delta$x),\n(pix)", fontsize = 15.0)
    ax4.axes.get_xaxis().set_visible(False)
    ax4.tick_params(labelsize=14)

    # y
    ax5 = fig.add_subplot(gs[5, :])
    ax5.axvline(x=400, linestyle=":", color="gray")
    ax5.axvline(x=1200, linestyle=":", color="gray")
    #ax5.plot(elapsed_time_inj, tilt_inj_simple_mod, color = "b", linewidth=10, alpha=0.15)
    ax5.plot(tilt_inj_positive_side, color = "b", linewidth=10, alpha=0.25)
    ax5.plot(tilt_inj_negative_side, color = "b", linewidth=10, alpha=0.25)
    ax5.plot(tilt_ret, color="red", marker='o', linestyle='None', markersize=2, alpha=1)
    ax5.axhline(y=10.7, linestyle = "--")
    ax5.axhline(y=-10.7, linestyle = "--")
    ax5.annotate("$+PS$", xy=(0, 0), xytext=(40, 11.5), textcoords="data", fontsize = 15.0)
    ax5.annotate("$-PS$", xy=(0, 0), xytext=(40, -13), textcoords="data", fontsize = 15.0)
    #ax5.set_title("Tilt, wrapped (x)", fontsize = 18.0)
    ax5.set_xlim([0,10000])
    ax5.set_ylim([-14,14])
    ax5.set_ylabel("Translation ($\Delta$y),\n(pix)", fontsize = 15.0)
    ax5.set_xlabel("Consecutive frame num (unitless)", fontsize = 15.0)
    ax5.tick_params(labelsize=14)

    plt.suptitle(csv_file)

    #plt.rcParams.update({'font.size': 12})

    #plt.tight_layout()

    plt.savefig("test.pdf")
    #plt.show()
