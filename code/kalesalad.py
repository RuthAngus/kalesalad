# Uses acf method to measure rotation periods for downloaded everest light
# curves.

import numpy as np
import matplotlib.pyplot as plt
import pyfits
from Kepler_ACF import corr_run
import os
from simple_acf import simple_acf
import sys
from multiprocessing import Pool
import pandas as pd
import glob
import scipy.stats as sps
import rotation as ro
import datetime

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 20,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def process_data(file):
    """
    Read the lightcurve from the fits format and sigma clip.
    prefix (str): the 4 digit number at the beginning of the epic id, e.g.
        "2011".
    id (str): the 4 digit number at the end of the epic id, e.g. "26368".
    c (str): campaign. e.g. "01"
    """
    with pyfits.open(file) as hdulist:
        time, flux = hdulist[1].data["TIME"], hdulist[1].data["FLUX"]
        # out = hdulist[1].data["OUTLIER"]

    m = np.isfinite(time) * np.isfinite(flux) #* (out < 1)
    x, med = time[m], np.median(flux[m])
    y = flux[m]/med - 1  # median normalise
    yerr = np.ones_like(y) * 1e-5

    # Sigma clip
    filtered_y = sps.sigmaclip(y, 2)[0]

    m = np.nonzero(np.in1d(y, filtered_y))[0]
    return x[m], y[m], yerr[m]


def run_acf(c, epic, plot=True):
    """
    Run the ACF on a light curve in the specified campaign.
    FOR PARALLEL RUNS.
    c (str): campaign, e.g. "c01".
    fn (str): fits file name for a target in campaign c.
    """

    #period, acf_smooth, lags, rvar, peaks, dips, leftdips, rightdips, \
    #bigpeaks = simple_acf(x, y)

    v = "2.0"
    filen = "hlsp_everest_k2_llc_{0}-c{1}_kepler_v{2}_lc.fits"\
        .format(epic, c.zfill(2), v)
    file = "data/c{0}/{1}".format(c.zfill(2), filen)

    # Load time and flux
    if not os.path.exists(file):
        print(file, "file not found")
        return None
    try:
        x, y, yerr = process_data(file)
    except (IOError, ValueError):
        print("Bad file", file)
        return None

    # compute the acf
    period, acf_smooth, lags, rvar, peaks = simple_acf(x, y)

    # make a plot
    if plot:
        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(x-x[0], y, "k.")
        plt.xlim(0, max(lags))
        plt.xlabel("$\mathrm{Time~(days)}$")
        plt.ylabel("$\mathrm{Normalised~flux}$")
        plt.subplot(2, 1, 2)
        plt.plot(lags, acf_smooth, "k")
        plt.xlabel("$\mathrm{lags~(days)}$")
        plt.ylabel("$\mathrm{ACF}$")
        plt.axvline(period, color="m")
        plt.xlim(min(lags), max(lags))
        plt.subplots_adjust(left=.16, bottom=.12, hspace=.4)
        plt.savefig("acfs/{}_acf".format(epic))

    # Measure LS period
    star = ro.prot(kepid=epic, x=x, y=y, yerr=yerr)
    pgram_period = star.pgram_ps(plot=True)

    return epic, period


def run_kalesalad(N):
    """
    Measure all rotation periods in a campaign - non parallel (for tests).
    """
    todays_date = datetime.date.today()
    results_file = "c{0}_periods_{1}.txt".format(c, todays_date)
    f.write("{0} {1} {2} {3}\n".format("epic_id", "ACF_period",
                                       "pgram_period", "pgram_period_err")

    c = str(sys.argv[1])
    df = pd.read_csv("c{}_targets.txt".format(c.zfill(2)), dtype=str)
    fns = df["epid"].values
    acf_periods, pgram_periods, pgram_period_errs, epics = [np.zeros(N) for i
                                                            in range(4)]
    for i, epic in enumerate(fns[:N]):
        v = "2.0"
        filen = "hlsp_everest_k2_llc_{0}-c{1}_kepler_v{2}_lc.fits"\
            .format(epic, c.zfill(2), v)
        file = "data/c{0}/{1}".format(c.zfill(2), filen)

        # Load time and flux
        if not os.path.exists(file):
            print(file, "file not found")
            return None
        try:
            x, y, yerr = process_data(file)
        except (IOError, ValueError):
            print("Bad file", file)
            return None

        # Measure ACF period
        _, acf_period = run_acf(c, epic, plot=True)

        # Measure LS period
        star = ro.prot(kepid=epic, x=x, y=y, yerr=yerr)
        pgram_period = star.pgram_ps(plot=True)

        with open(results_file, "a") as f:
            f.write("{0} {1} {2} {3}\n".format(epic, acf_period,
                                               pgram_period[0],
                                               pgram_period[1]))

if __name__ == "__main__":
    from functools import partial
    c = str(sys.argv[1])

    # open("c{0}_periods.txt".format(c), "w")

    run_kalesalad(100)
    assert 0

    df = pd.read_csv("c{}_targets.txt".format(c.zfill(2)), dtype=str)
    fns = df["epid"].values

    f = partial(run_acf, c)

    pool = Pool()
    for val in pool.map(f, fns):
        if val is None:
            continue
        epic, acf_period, epic_period = val
        # append data to file
        with open("c{0}_periods.txt".format(c), "a") as f:
            f.write("{0} {1} \n".format(epic, period))
