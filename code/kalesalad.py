from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from Kepler_ACF import corr_run
import os
from simple_acf import simple_acf
import sys
from multiprocessing import Pool

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

    hdulist = pyfits.open(file)
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

    if os.path.exists(file):
        x, y = process_data(file)

        # compute the acf
        period, acf_smooth, lags, rvar, peaks, dips, leftdips, rightdips, \
                bigpeaks = simple_acf(x, y)

        # append data to file
        with open("c{0}_periods.txt".format(c), "a") as f:
            f.write("{0} {1} \n".format(epic, period))

        # make a plot
        if plot:
            plt.clf()
            plt.subplot(2, 1, 1)
            plt.plot(x, y, "k.")
            plt.xlim(min(x), max(x))
            plt.xlabel("$\mathrm{Time~(days)}$")
            plt.ylabel("$\mathrm{Normalised~flux}$")
            plt.subplot(2, 1, 2)
            plt.plot(lags, acf_smooth, "k")
            plt.xlabel("$\mathrm{lags~(days)}$")
            plt.ylabel("$\mathrm{ACF}$")
            plt.axvline(p, color="m")
            plt.savefig("results/{}_acf".format(epic))
    else:
        print(file, "file not found")


def run_kalesalad_multi(index):
    """
    Measure all rotation periods in a campaign using parallel processing.
    Iterate over fits file names.
    """
    c = str(sys.argv[1])
    fns = np.genfromtxt("c{0}_targets.txt".format(c), dtype=str).T
    run_acf(c, fns[index], plot=False)

def run_kalesalad(N):
    """
    Measure all rotation periods in a campaign - non parallel (for tests).
    """
    c = str(sys.argv[1])
    fns = np.genfromtxt("c{0}_targets.txt".format(c), dtype=str).T
    for fn in fns[:N]:
        run_acf(c, fn, plot=False)


if __name__ == "__main__":
    c = str(sys.argv[1])

    assert os.path.exists("c{0}_periods.txt".format(c)) == False, \
            "You need to delete the old file!"

    fns = np.genfromtxt("c{0}_targets.txt".format(c), dtype=str).T
    run_kalesalad(10)

#     pool = Pool()
#     pool.map(run_kalesalad_multi, range(len(fns)))
