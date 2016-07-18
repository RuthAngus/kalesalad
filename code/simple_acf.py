import numpy as np
import matplotlib.pyplot as plt
import emcee
import glob

def simple_acf(x, y):

    # fit and subtract straight line
    AT = np.vstack((x, np.ones_like(x)))
    ATA = np.dot(AT, AT.T)
    m, b = np.linalg.solve(ATA, np.dot(AT, y))
    y -= m*x + b

    # perform acf
    acf = dan_acf(y)

    # create 'lags' array
    gap_days = 0.02043365
    lags = np.arange(len(acf))*gap_days

    N = len(acf)
    double_acf, double_lags = [np.zeros((2*N)) for i in range(2)]
    double_acf[:N], double_lags[:N] = acf[::-1], -lags[::-1]
    double_acf[N:], double_lags[N:] = acf, lags
    acf, lags = double_acf, double_lags

    # smooth with Gaussian kernel convolution
    Gaussian = lambda x, sig: 1./(2*np.pi*sig**.5) * np.exp(-0.5*(x**2)/
                                                            (sig**2))
    conv_func = Gaussian(np.arange(-28, 28, 1.), 9.)
    acf_smooth = np.convolve(acf, conv_func, mode='same')
    acf_smooth, lags = acf_smooth[N:], lags[N:]

    # ditch the first point
    acf_smooth, lags = acf_smooth[1:], lags[1:]

    peaks, dips, leftdips, rightdips, bigpeaks = find_peaks(acf_smooth, lags)

    # find the first and second peaks
    if len(peaks) > 1:
        if acf_smooth[peaks[0]] > acf_smooth[peaks[1]]:
            period = lags[peaks[0]]
        else:
            period = lags[peaks[1]]
    elif len(peaks) == 1:
        period = lags[peaks][0]
    elif not len(peaks):
        period = np.nan

    # find the highest peak
    m = acf_smooth == max(acf_smooth[peaks])
    highest_peak = acf_smooth[m][0]
    period = lags[m][0]
    print(highest_peak, period)

    rvar = np.percentile(y, 95)

    return period, acf_smooth, lags, rvar, peaks, dips, leftdips, rightdips, \
        bigpeaks


def find_peaks(acf_smooth, lags, t=.2):

    # find all the peaks
    peaks = np.array([i for i in range(1, len(lags)-1)
                     if acf_smooth[i-1] < acf_smooth[i] and
                     acf_smooth[i+1] < acf_smooth[i]])

    # find all the dips
    dips = np.array([i for i in range(1, len(lags)-1)
                     if acf_smooth[i-1] > acf_smooth[i] and
                     acf_smooth[i+1] > acf_smooth[i]])

    # find all dips immediately left and right of the peaks
    leftdips = np.array([lags[dips][lags[dips] < lags[peak]][-1] for peak in
                         peaks])
    if lags[dips[-1]] > lags[peaks[-1]]:
        rightdips = np.array([lags[dips][lags[dips] > lags[peak]][0] for peak
                              in peaks])
    else:  # if there is no right hand dip:
        rightdips = [lags[dips][lags[dips] > lags[peaks[i]]][0]
                     for i in range(len(peaks)-1)]
        rightdips.append(lags[-1])
        rightdips = np.array(rightdips)

    # find the peak heights
    leftdip_heights = acf_smooth[lags == leftdips]
    rightdip_heights = acf_smooth[lags == rightdips]
    peak_heights = acf_smooth[peaks]
    meandiffs = .5*(np.abs(peak_heights - leftdip_heights) +
                    np.abs(peak_heights - rightdip_heights))
    bigpeaks = peaks[meandiffs > t]

    return peaks, dips, leftdips, rightdips, bigpeaks


def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


# dan's acf function
def dan_acf(x, axis=0, fast=False):
    """
    Estimate the autocorrelation function of a time series using the FFT.
    :param x:
        The time series. If multidimensional, set the time axis using the
        ``axis`` keyword argument and the function will be computed for every
        other axis.
    :param axis: (optional)
        The time axis of ``x``. Assumed to be the first axis if not specified.
    :param fast: (optional)
        If ``True``, only use the largest ``2^n`` entries for efficiency.
        (default: False)
    """
    x = np.atleast_1d(x)
    m = [slice(None), ] * len(x.shape)

    # For computational efficiency, crop the chain to the largest power of
    # two if requested.
    if fast:
        n = int(2**np.floor(np.log2(x.shape[axis])))
        m[axis] = slice(0, n)
        x = x
    else:
        n = x.shape[axis]

    # Compute the FFT and then (from that) the auto-correlation function.
    f = np.fft.fft(x-np.mean(x, axis=axis), n=2*n, axis=axis)
    m[axis] = slice(0, n)
    acf = np.fft.ifft(f * np.conjugate(f), axis=axis)[m].real
    m[axis] = 0
    return acf / acf[m]


if __name__ == "__main__":

    DIR = "."  # edit me!
    fnames = glob.glob("%s/*.dat" % DIR)

    for i, fname in enumerate(fnames[1:]):
        id = fname.split("/")[-1].split("_")[0]  # edit me!
        x, y, _, _ = np.genfromtxt(fname, skip_header=1).T
        yerr = np.ones_like(y) * 1e-5  # FIXME

        period, acf, lags = simple_acf(x, y)
        make_plot(acf, lags, id)
