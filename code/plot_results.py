from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from download_epic import get_catalog
import os
import sys
from matplotlib.ticker import NullFormatter


def period2age(period, bv):
    a, b, c, n = .4, .31, .4, .55
    return (period / (a*(bv-c)**b) )**(1./n) * 1e-3


def get_EB_catalog(name, basepath="data"):
    fn = os.path.join(basepath, "{0}.h5".format(name))
    if os.path.exists(fn):
        return pd.read_hdf(fn, name)
    data = np.genfromtxt(os.path.join(basepath, "EB_catalog_periods.txt"),
                         skip_header=65).T
    EBids, EBperiods = data[0], data[3]
    with open(os.path.join(basepath, "{0}.csv".format(name)), "a") as f:
        for i in range(len(EBids)):
            f.write("{0}, {1}".format(EBids[i], EBperiods[i]))
    csvfn = os.path.join(basepath, "{0}.csv".format(name))
    df = pd.read_csv(csvfn, low_memory=False)
    df.to_hdf(fn, name, format="t")
    return df


def match(myids, myperiods, df, no_binaries=True):
    """
    Match the rotation periods with information from the epic.
    Return a dictionary with matched information
    """

    bv, bverr, teff, tefferr1, tefferr2, ra, raerr, dec, decerr, logg, \
            loggerr1, loggerr2, lon, lonerr, lat, laterr = \
            [np.zeros_like(myids) for i in range(16)]

    if no_binaries:
        ebids = remove_binaries()  # load the list of binary ids

    for i in range(len(myids)):

        if no_binaries:  # if the target is a binary, skip it (leave a zero)
            if any(myids[i] == ebids):
                bv[i], bverr[i], teff[i], tefferr1[i] = 0, 0, 0, 0
                tefferr2[i], ra[i], raerr[i], dec[i], decerr[i] = 0, 0, 0, 0, 0
                logg[i], loggerr1[i], loggerr2[i] = 0, 0, 0

        m = myids[i] == df.epic_number
        bv[i] = df.k2_bjmag[m] - df.k2_vjmag[m]
        bverr[i] = (df.k2_bjmagerr[m]**2 + df.k2_vjmagerr[m]**2)**.5
        teff[i] = df.k2_teff[m]
        tefferr1[i], tefferr2[i] = df.k2_tefferr1[m], df.k2_tefferr2[m]
        ra[i], raerr[i] = df.ra[m], df.raerr[m]
        dec[i], decerr[i] = df.dec[m], df.decerr[m]
        logg[i] = df.k2_logg[m]
        loggerr1[i], loggerr2[i] = df.k2_loggerr1[m], df.k2_loggerr2[m]
        lon[i], lonerr[i] = df.k2_glon[m], df.k2_glonerr[m]
        lat[i], laterr[i] = df.k2_glat[m], df.k2_glaterr[m]

    return {"id": myids,
            "period": myperiods,
            "bv": bv,
            "bverr": bverr,
            "teff": teff,
            "tefferr1": tefferr1,
            "tefferr2": tefferr2,
            "ra": ra,
            "raerr": raerr,
            "dec": dec,
            "decerr": decerr,
            "logg": logg,
            "loggerr1": loggerr1,
            "loggerr2": loggerr2,
            "lon": lon,
            "lonerr": lonerr,
            "lat": lat,
            "laterr": laterr}


def remove_binaries(path="data"):
    data = np.genfromtxt(os.path.join(path, "EB_catalog_class.txt"),
                         skip_header=58).T
    ids, detached_prob, non_detached_prob = data[0], data[4], data[5]
    m = (detached_prob > np.median(detached_prob)) * \
            (non_detached_prob > np.median(non_detached_prob))
    return ids[m]  # array of binary ids

def combine_campaign_data(c, prefixes):
    myids, periods, err = [], [], []
    for p in prefixes:
        data = np.genfromtxt("c{0}_{1}_periods.txt".format(c, p)).T
        myids.append(data[0])
        periods.append(data[1])
        err.append(data[2])
    return np.array([i for j in myids for i in j]), \
            np.array([i for j in periods for i in j]), \
            np.array([i for j in err for i in j])


if __name__ == "__main__":

    c = str(sys.argv[1])

#     if c == "01":  # combine all campaign 1 data
#         prefixes = ["2011", "2012", "2013", "2014", "2015", "2016", "2017",
#                     "2018", "2019", "2102"]
#         myids, periods, err = combine_campaign_data(c, prefixes)
#     else:

    myids, periods = np.genfromtxt("c{0}_periods.txt".format(c)).T

    # plot armstrong data
#     data = np.genfromtxt("data/EB_catalog_periods.txt", skip_header=65).T
#     EBids, EBperiods = data[0], data[3]
#     k2 = match(EBids, EBperiods, df)

    df = get_catalog("k2targets")
#     print(df.keys())

    k2 = match(myids, periods, df)

    plotpar = {'axes.labelsize': 20,
               'text.fontsize': 20,
               'legend.fontsize': 20,
               'xtick.labelsize': 20,
               'ytick.labelsize': 20,
               'text.usetex': True}
    plt.rcParams.update(plotpar)

#     plt.clf()
#     plt.plot(k2["teff"], k2["logg"], "k.")
#     print("plotting", len(k2["teff"]), "teffs")
#     plt.xlabel("$\mathrm{T}_{\mathrm{eff}}~\mathrm{(K)}$")
#     plt.ylabel("$\log~(g)~[\log_{10}\mathrm{(cgs)}]$")
#     plt.xlim(7500, 1800)
#     plt.ylim(6, 0)
#     plt.savefig("figs/teff_vs_logg_c{0}".format(c))
#
#     plt.clf()
#     age = period2age(k2["period"], k2["bv"])
# #     m = (k2["period"] < 45) * (k2["period"] > 0)  # take out bad periods
#     # take out hot stars, bad periods and bad ages.
#     m = (k2["period"] < 45) * (k2["period"] > 0)
# #     * (k2["bv"]>.4) * (age > 0) \
# #             * (age < 13) * np.isfinite(age) * (k2["logg"] > 4.2)
#     plt.scatter(k2["ra"][m], k2["dec"][m], marker="o", c=k2["period"][m],
#             edgecolor="", cmap="viridis", s=5, vmin=0, vmax=40)
#     print("plotting", len(k2["dec"][m]), "decs")
#     plt.ylabel("$\mathrm{Declination~(degrees)}$")
#     plt.xlabel("$\mathrm{Right~Ascension~(degrees)}$")
#     plt.colorbar(label="$\mathrm{P}_{\mathrm{rot}}~\mathrm{(days)}$")
#     plt.subplots_adjust(bottom=.2)
#     middec = .5*(min(k2["dec"]) + max(k2["dec"]))
#     plt.ylim(middec - 9, middec + 9)
#     plt.savefig("figs/ra_vs_dec_period_c{0}".format(c))
#
#     plt.clf()
#     plt.hist(k2["period"][m], histtype="stepfilled", color="w")
#     plt.ylabel("$\mathrm{N}$")
#     plt.xlabel("$\mathrm{P}_{\mathrm{rot}}~\mathrm{(days)}$")
# #     plt.title("Lat = {0}, long = {1}".format(.5*(min(k2["lat"]) +
# #               max(k2["lat"])), .5*(min(k2["lon"]) + max(k2["lon"]))))
#     plt.xlim(0, 40)
#     plt.savefig("figs/period_hist_c{0}".format(c))
#
#     plt.clf()
#     # take out hot stars, bad periods and bad ages.
#     m = (k2["period"] < 45) * (k2["period"] > 0) * (k2["bv"]>.4) * (age > 0) \
#             * (age < 13) * np.isfinite(age) * (k2["logg"] > 4.2)
#     plt.scatter(k2["dec"][m], k2["ra"][m], marker="o", c=age[m],
#             edgecolor="", cmap="viridis", s=3)
#     print("plotting", len(k2["dec"][m]), "decs")
#     plt.xlabel("$\mathrm{Declination~(degrees)}$")
#     plt.ylabel("$\mathrm{Right~Ascension~(degrees)}$")
#     plt.colorbar(label="$\mathrm{Age~(Gyr)}$")
#     plt.subplots_adjust(bottom=.2)
#     plt.xlim(middec - 9, middec + 9)
#     plt.savefig("figs/ra_vs_dec_age_c{0}".format(c))
#
#     plt.clf()
#     plt.hist(age[m], histtype="stepfilled", color="w")
#     plt.ylabel("$\mathrm{N}$")
#     plt.xlabel("$\mathrm{Age (Gyr)}$")
#     plt.xlim(0, 13.7)
#     plt.title("Lat = {0}, long = {1}".format(.5*(min(k2["lat"]) +
#               max(k2["lat"])), .5*(min(k2["lon"]) + max(k2["lon"]))))
#     plt.savefig("figs/age_hist_c{0}".format(c))
#
#
#     plt.clf()
#     # take out hot stars, bad periods and bad ages.
#     m = (k2["period"] < 45) * (k2["period"] > 0) * (k2["bv"]>.4) * (age > 0) \
#             * (age < 13) * np.isfinite(age) * (k2["logg"] > 4.2)
#     plt.scatter(k2["lon"][m], k2["lat"][m], marker="o", c=age[m],
#                 edgecolor="", cmap="RdPu", s=3)
#     print("plotting", len(k2["lat"][m]), "lats")
#     plt.xlabel("$\mathrm{Galactic~longitude}$")
#     plt.ylabel("$\mathrm{Galactic~latitude}$")
#     plt.colorbar(label="$\mathrm{Age~(Gyr)}$")
#     plt.subplots_adjust(bottom=.2)
#     plt.savefig("figs/lon_vs_lat_c{0}".format(c))

    plt.clf()
    left, width = .1, .65
    bottom, height = .1, .65
    bottom_h = left_h = left + width + .02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, .2]
    rect_histy = [left_h, bottom, .2, height]
    plt.figure(1, figsize=(8, 8))
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # take out bad periods and giants.
    m = (k2["period"] < 45) * (k2["period"] > 0) * (k2["logg"] > 4.2)
    p = k2["period"][m]
    teff = k2["teff"][m]
    axScatter.scatter(teff, p, c="k", s=3)
    axScatter.set_xlabel("$\mathrm{T}_{\mathrm{eff}}~\mathrm{(K)}$")
    axScatter.set_ylabel("$\mathrm{P}_{\mathrm{rot}}~\mathrm{(days)}$")
    axScatter.set_xlim(2000, 7000)
    axScatter.set_ylim(0, 45)
    print("plotting", len(p), "periods")

#     binwidth = .25
#     xymax = np.max( [np.max(np.fabs(teff)), np.max(np.fabs(p))] )
#     lim = (int(xymax/binwidth) + 1) * binwidth
#     axScatter.set_xlim((-lim, lim))
#     axScatter.set_ylim((-lim, lim))
#     bins = np.arange(-lim, lim + binwidth, binwidth)
#     axHistx.hist(teff, bins=bins)
#     axHisty.hist(p, bins=bins, orientation="horizontal")
#     axHistx.set_xlim(axScatter.get_xlim())
#     axHisty.set_ylim(axScatter.get_ylim())
    axHistx.hist(teff, histtype="stepfilled", color="w")
    axHisty.hist(p, histtype="stepfilled", color="w", orientation="horizontal")

#     plt.xlim(max(teff), min(teff))
#     plt.xlabel("$\mathrm{T}_{\mathrm{eff}}~\mathrm{(K)}$")
#     plt.ylabel("$\mathrm{P}_{\mathrm{rot}}~\mathrm{(days)}$")
#     plt.yscale("log")
#     plt.subplots_adjust(left=.2)
    plt.savefig("figs/period_vs_bv_c{0}".format(c))
