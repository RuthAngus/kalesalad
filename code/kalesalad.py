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

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 20,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def process_data(file):
    """
    Read the lightcurve from the fits format.
    prefix (str): the 4 digit number at the beginning of the epic id, e.g.
        "2011".
    id (str): the 4 digit number at the end of the epic id, e.g. "26368".
    c (str): campaign. e.g. "01"
    """
    with pyfits.open(file) as hdulist:
        time, flux = hdulist[1].data["TIME"], hdulist[1].data["FLUX"]
        out = hdulist[1].data["OUTLIER"]

    m = np.isfinite(time) * np.isfinite(flux) * (out < 1)
    x, med = time[m], np.median(flux[m])
    y = flux[m]/med - 1  # median normalise
    return x, y


def run_acf(c, fn, plot=False):
    """
    Run the ACF on a light curve in the specified campaign.
    c (str): campaign, e.g. "c01".
    fn (str): fits file name for a target in campaign c.
    """

    epic = fn[20:29]
    file = "data/c{0}/{1}".format(c, fn)

    if not os.path.exists(file):
        print(file, "file not found")
        return None

    try:
        x, y = process_data(file)

#             period, acf_smooth, lags, rvar, peaks, dips, leftdips, rightdips, \
#                     bigpeaks = simple_acf(x, y)

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
        plt.savefig("results/{}_acf".format(epic))

    return epic, period


def run_kalesalad(N):
    """
    Measure all rotation periods in a campaign - non parallel (for tests).
    """
    c = str(sys.argv[1])
    df = pd.read_csv("c{}_targets.txt".format(c.zfill(2)), dtype=str)
    fns = df["epid"].values
    for fn in fns[:N]:
        run_acf(c, fn, plot=False)


if __name__ == "__main__":
    from functools import partial
    c = str(sys.argv[1])

    open("c{0}_periods.txt".format(c), "w")

    # fns = np.genfromtxt("c{0}_targets.txt".format(c), dtype=str).T
    print("c{}_targets.txt".format(c.zfill(2)))

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
