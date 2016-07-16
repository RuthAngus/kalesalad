import numpy as np
import matplotlib.pyplot as plt
from simple_acf import simple_acf
import pyfits
from Kepler_ACF import corr_run
import tarfile
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
    # import everest
    # everest.compute.Compute(201126368)
    # assert 0

    # load lightcurves
    ids = np.genfromtxt("ids.txt", dtype=str)
    prefixes = np.genfromtxt("prefixes.txt", dtype=str)
    c = "01"
    for p in prefixes:
        for id in ids:

            prefix = p[:4]
            epic = "{0}{1}".format(prefix, id)
            path = "data/c01/{0}00000/{1}".format(prefix, id)
            end = "kepler_v1.0_lc.fits"
            file = "{0}/hlsp_everest_k2_llc_{1}-c{2}_{3}".format(path, epic,
                                                                 c, end)
            tf = "data/c{0}_{1}.tar.gz".format(c, p)
            with tarfile.open(tf) as f:
                print(f)
                files = f.getmember(file)
                print(files)
            assert 0

            if os.path.exists(file):
                x, y = process_data(file)
                plt.clf()
                plt.plot(x, y, "k.")
                plt.savefig("results/{}_lc".format(epic))

                period, acf_smooth, lags, rvar = simple_acf(x, y)
                plt.clf()
                plt.plot(lags, acf_smooth)
                plt.axvline(period, color="r")
                plt.title("P = {0:.2f} Rvar = {1:.2f}".format(period, rvar))
                plt.savefig("results/{}_acf".format(epic))
