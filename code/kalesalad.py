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


def run_acf(c, p, ids, plot=False):
    """
    Run the ACF on all light curves in the specified campaign, with the
    specified prefix.
    c (str): campaign, e.g. "c01".
    p (int): prefix, e.g. 2011.
    ids (np.array of ints): list of epic ids (9 digits).
    """

    assert os.path.exists("c{0}_{1}_periods.txt".format(c, str(p)[:4])) == \
            False, "You need to delete the old file!"

    for i, id in enumerate(ids):

        prefix = str(p)[:4]
        epic = "{0}{1}".format(prefix, id)
        print(epic)
        path = "data/c01/{0}00000/{1}".format(prefix, id)
        end = "kepler_v1.0_lc.fits"
        file = "{0}/hlsp_everest_k2_llc_{1}-c{2}_{3}".format(path, epic,
                                                                c, end)

        if os.path.exists(file):
            print(i, "of", len(ids))
            x, y = process_data(file)

            # compute the acf
#             acf_smooth, lags, period, err, locheight = corr_run(x, y)
            period, acf_smooth, lags, rvar, peaks, dips, leftdips, rightdips, \
                    bigpeaks = simple_acf(x, y)
            err = 0

            # append data to file
            with open("c{0}_{1}_periods.txt".format(c, prefix), "a") as f:
                f.write("{0} {1} {2} \n".format(epic, period, err))

            # make a plot
            if plot:
                p = period
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
                plt.axvline(p, color="CornFlowerBlue",
                            label="${0:.2f}$".format(p))
                plt.legend(loc="best")
                plt.savefig("results/{}_acf".format(epic))


def run_kalesalad(index):
    """
    Measure all rotation periods in campaign 1 using parallel processing
    Campaign 1:
    p = ["2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019",
           "2102"]
    c = "01"
    """
    # Campaign 1
    p = ["2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019",
           "2102"]
    p = ["2011"]
    c = "01"
    # Campaign 2
#     p = ["2024", "2025", "2026", "2027", "2028", "2029", "2030", "2031",
#            "2032", "2033", "2034", "2035", "2036", "2037", "2038", "2039",
#            "2040", "2041", "2042", "2043", "2044", "2045", "2046", "2047",
#            "2048", "2049", "2050", "2051", "2052", "2053", "2054", "2056",
#            "2057"]
#     c = "02"
    ids = np.genfromtxt("{0}_ids.txt".format(p[index]), dtype=str)
    run_acf(c, p[index], ids, plot=True)


if __name__ == "__main__":
    pool = Pool()
    pool.map(run_kalesalad, range(1))
