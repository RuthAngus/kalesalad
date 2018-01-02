# Find the giants in the TGAS/EPIC crossmatched catalog and remove them from
# the list.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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


if __name__ == "__main__":

    # Plot a CMD.
    df = pd.read_csv("epic_tgas.csv")

    # Calculate colours and magnitudes
    abs_jmag = abs_mag(df.k2_jmag.values, df.tgas_parallax.values)
    colour = df.k2_jmag.values - df.k2_kmag.values

    # Remove NaNs
    m = np.isfinite(colour) * np.isfinite(abs_jmag)
    colour, abs_jmag = colour[m], abs_jmag[m]

    counts, xbins, ybins = np.histogram2d(colour, abs_jmag, bins=100)

    plt.clf()
    plt.scatter(colour, abs_jmag, c=colour, s=2)
    # plt.contour(xbins[:-1], ybins[:-1], counts, colors='black');
    plt.contour(counts.transpose(),
                extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                linewidths=.5, colors='white', linestyles='solid')
    plt.colorbar(label="$J - K$")
    plt.xlabel("$J - K$")
    plt.ylabel("$M_J$")
    plt.ylim(22, 8)
    plt.xlim(-.2, 1.5)
    plt.subplots_adjust(left=.15, bottom=.15)
    plt.savefig("CMD")
