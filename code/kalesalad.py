# Uses acf method to measure rotation periods for downloaded everest light
# curves.

from __future__ import print_function
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


def run_acf(c, epic, x, y, plot=False):
    """
    Run the ACF on a light curve in the specified campaign.
    c (str): campaign, e.g. "c01".
    fn (str): fits file name for a target in campaign c.
    """

    #period, acf_smooth, lags, rvar, peaks, dips, leftdips, rightdips, \
    #bigpeaks = simple_acf(x, y)

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

    return epic, period


def run_kalesalad(N):
    """
    Measure all rotation periods in a campaign - non parallel (for tests).
    """
    c = str(sys.argv[1])
    df = pd.read_csv("c{}_targets.txt".format(c.zfill(2)), dtype=str)
    fns = df["epid"].values
    # periods, epics = [np.zeros(N) for i in range(2)]
    for i, epic in enumerate(fns[:N]):

        v = "2.0"
        filen = "hlsp_everest_k2_llc_{0}-c{1}_kepler_v{2}_lc.fits"\
            .format(epic, c.zfill(2), v)
        print(filen)
        file = "data/c{0}/{1}".format(c.zfill(2), filen)
        print(file)

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
        print(epic)
        _, acf_period = run_acf(c, epic, x, y, plot=True)

        # Measure LS period
        star = ro.prot(kepid=epic, x=x, y=y, yerr=yerr)
        pgram_period = star.pgram_ps(plot=True)
        print(pgram_period)

        # periods[i] = period
        # epics[i] = epic
    # res = dict({"epic": np.array(epics), "period": np.array(periods)})
    # res = pd.DataFrame(res)
    # res.to_csv("results.")

if __name__ == "__main__":
    from functools import partial
    c = str(sys.argv[1])

    open("c{0}_periods.txt".format(c), "w")

    run_kalesalad(2)
    assert 0

    f = partial(run_acf, c)

    pool = Pool()
    for val in pool.map(f, fns):
        if val is None:
            continue
        epic, period = val
        # append data to file
        with open("c{0}_periods.txt".format(c), "a") as f:
            f.write("{0} {1} \n".format(epic, period))
