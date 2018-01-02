# Find the giants in the TGAS/EPIC crossmatched catalog and remove them from
# the list.
# Saves a .csv file of the tgas_epic.csv catalogue with the giants removed and
# stars of very high and low temperatures removed.
# I have not made a cut in parallax error.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import isochrones
# from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones.mist import MIST_Isochrone
from scipy.interpolate import interp1d

plotpar = {'axes.labelsize': 20,
            'text.fontsize': 20,
            'legend.fontsize': 20,
            'xtick.labelsize': 20,
            'ytick.labelsize': 20,
            'text.usetex': True}
plt.rcParams.update(plotpar)

def abs_mag_w_MC_uncertainty(mean_m, m_err, mean_p, parallax_err, N):
    values = np.vstack((mean_m + np.random.randn(N)*m_err,
                        mean_p + np.random.randn(N)*parallax_err)).T
    abs_mags = [abs_mag(m, parallax) for m, parallax in values]
    mean = np.mean(abs_mags)
    lower, upper = np.percentile(abs_mags, 16), np.percentile(abs_mags, 84)
    return mean, mean - lower, upper - mean


def abs_mag(m, parallax):
    return m - 5 * np.log10(1./parallax) + 5

def CMD_plot(colour, abs_mag):
    """
    Plot and save a CMD.
    """

    counts, xbins, ybins = np.histogram2d(colour, abs_jmag, bins=100)

    plt.clf()
    lower_lim, upper_lim = -.2, 1.5
    m = (lower_lim < colour) * (colour < upper_lim)
    plt.scatter(colour[m], abs_jmag[m], c=colour[m], s=2)
    # plt.contour(xbins[:-1], ybins[:-1], counts, colors='black');
    plt.contour(counts.transpose(),
                extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                linewidths=.5, colors='white', linestyles='solid')
    plt.colorbar(label="$J - K$")
    plt.xlabel("$J - K$")
    plt.ylabel("$M_J$")
    plt.ylim(22, 8)
    plt.xlim(lower_lim, upper_lim)
    plt.subplots_adjust(left=.15, bottom=.15)
    plt.savefig("CMD")


def HRD_plot(teff, abs_gmag):

    # dar = Dartmouth_Isochrone()
    mist = MIST_Isochrone()
    iso_300 = mist.isochrone(age=np.log10(300e6), feh=0.0, AV=0.0)

    counts, xbins, ybins = np.histogram2d(teff, abs_gmag, bins=100)

    # Make temperature cuts
    lower_lim, upper_lim = 3200, 8000
    m = (lower_lim < teff) * (teff < upper_lim)

    plt.clf()
    # Plot points + contours
    plt.scatter(teff[m], abs_gmag[m], c=teff[m], s=2, cmap="viridis_r",
                zorder=0)
    plt.colorbar(label="$T_{\mathrm{eff}}$")
    plt.contour(counts.transpose(),
                extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                linewidths=.5, colors='white', linestyles='solid')

    # Plot isochrones
    plt.plot(iso_300.Teff[:100], iso_300.G_mag[:100]+1, "k", lw=.5)
    plt.plot(iso_300.Teff[:100], iso_300.G_mag[:100]-1, "k", lw=.5)

    # Select stars between isochrones
    fup = interp1d(iso_300.Teff[:100], iso_300.G_mag[:100] + 1)
    flo = interp1d(iso_300.Teff[:100], iso_300.G_mag[:100] - 1)
    inds = []
    for i, G in enumerate(abs_gmag[m]):
        upper_G_at_this_teff = fup(teff[m][i])
        lower_G_at_this_teff = flo(teff[m][i])
        if (lower_G_at_this_teff < G) and (G < upper_G_at_this_teff):
            inds.append(i)
    # plt.scatter(teff[m][inds], abs_gmag[m][inds], c="k", s=2, zorder=1)

    plt.xlabel("$T_{\mathrm{eff}}$")
    plt.ylabel("$M_G$")
    plt.ylim(10, -7)
    plt.xlim(upper_lim, lower_lim)
    plt.subplots_adjust(left=.15, bottom=.15)
    plt.savefig("HRD")
    plt.savefig("HRD.pdf")

    return m, inds


if __name__ == "__main__":
    # import isochrones.dartmouth
    # isochrones.dartmouth.download_grids()

    # Plot a CMD.
    df = pd.read_csv("epic_tgas.csv")

    # Calculate colours and magnitudes
    abs_jmag = abs_mag(df.k2_jmag.values, df.tgas_parallax.values*1e-3)
    colour = df.k2_jmag.values - df.k2_kmag.values

    # Remove NaNs
    m = np.isfinite(colour) * np.isfinite(abs_jmag)
    colour, abs_jmag = colour[m], abs_jmag[m]

    CMD_plot(colour, abs_jmag)

    teff = df.k2_teff.values
    abs_gmag = abs_mag(df.tgas_phot_g_mean_mag.values,
                       df.tgas_parallax.values*1e-3)
    m = np.isfinite(teff) * np.isfinite(abs_gmag)
    teff, abs_gmag = teff[m], abs_gmag[m]

    teff_cut, inds = HRD_plot(teff, abs_gmag)

    df_temp_cuts = df.iloc[teff_cut]
    df_dwarfs = df_temp_cuts.iloc[inds]
    df_dwarfs.to_csv("tgas_epic_dwarfs.csv")
