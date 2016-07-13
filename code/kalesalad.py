import numpy as np
import matplotlib.pyplot as plt
from simple_acf import simple_acf
import pyfits
from Kepler_ACF import corr_run
import tarfile


def process_data(prefix, id, c):
    """
    Read the lightcurve from the fits format.
    prefix (str): the 4 digit number at the beginning of the epic id, e.g.
        "2011".
    id (str): the 4 digit number at the end of the epic id, e.g. "26368".
    c (str): campaign. e.g. "01"
    """

    path = "data/c01/{0}00000/{1}".format(prefix, id)
    file = "{0}/hlsp_everest_k2_llc_{1}-c{2}_kepler_v1.0_lc.fits".format(path,
                                                                      epic, c)

    # with open(file, "rb") as tarf:
    #     with tarfile.open(fileobj=tarf) as f:
    #         for mem in f.getmembers():
    #             if not fnmatch.fnmatch(mem.name, "*.fits"):
    #                 continue
    #             fitsf = f.extractfile(mem)

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
    prefixes = np.genfromtxt("prefixes.txt", dtype=str)
    for p in prefixes:
        for id in ids:

            prefix = p[:4]
            epic = "{0}{1}".format(prefix, id)
            c = "01"
            x, y = process_data(prefix, id, c)
            plt.clf()
            plt.plot(x, y, "k.")

            period, acf_smooth, lags = simple_acf(x, y)
            plt.clf()
            plt.plot(lags, acf_smooth)
            plt.axvline(period, color="r")
            plt.savefig("acfs/{}".format(epic))
