import numpy as np
import matplotlib.pyplot as plt
import pyfits
from Kepler_ACF import corr_run
import os


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


if __name__ == "__main__":

    # load lightcurves
    ids = np.genfromtxt("ids.txt", dtype=str)
    p = np.genfromtxt("prefixes.txt", dtype=str)
    c = "01"
    periods, perr, epics = [], [], []
    for i, id in enumerate(ids):

        prefix = str(p)[:4]
        epic = "{0}{1}".format(prefix, id)
        path = "data/c01/{0}00000/{1}".format(prefix, id)
        end = "kepler_v1.0_lc.fits"
        file = "{0}/hlsp_everest_k2_llc_{1}-c{2}_{3}".format(path, epic,
                                                                c, end)
        if os.path.exists(file):
            print(i, "of", len(ids))
            x, y = process_data(file)
            plt.clf()

            acf_smooth, lags, period, err, locheight = corr_run(x, y)
            plt.clf()
            plt.subplot(2, 1, 1)
            plt.plot(x, y, "k.")
            plt.subplot(2, 1, 2)
            plt.plot(lags, acf_smooth)
            plt.axvline(period, color="r")
            plt.title("P = {0:.2f}".format(period))
            plt.savefig("results/{}_acf".format(epic))
            periods.append(period)
            perr.append(err)
            epics.append(epic)

    data = np.vstack((np.array(epics), np.array(periods), np.array(perr)))
    np.savetxt("{0}{1}_periods.txt", data.T)
