import numpy as np
import matplotlib.pyplot as plt
from kalesalad import process_data
from mklc import mklc
from simple_acf import simple_acf


def k2lc(epic):
    """
    load k2 light curve
    """
    prefix = epic[:4]
    id = epic[4:]
    c = "01"
    path = "data/c01/{0}00000/{1}".format(prefix, id)
    end = "kepler_v1.0_lc.fits"
    file = "{0}/hlsp_everest_k2_llc_{1}-c{2}_{3}".format(path, epic, c, end)
    x, y = process_data(file)
    return x, y


def injection_params(N, params):
    """
    Generate parameter distributions.
    nspots - log uniform.
    inclination - uniform in sin^2(i)j.
    rotation period - log uniform.
    amplitude - log uniform.
    tau - (multiple of rotation period) log uniform
    """

    nspot_min, nspot_max, incl_min, incl_max, amp_min, amp_max, pmin, pmax, \
        tau_min, tau_max = params
    nspots = np.exp(np.random.uniform(np.log(nspot_min), np.log(nspot_max),
                                      N))
    incl = np.arcsin(np.random.uniform(np.sin(incl_min)**2,
                                       np.sin(incl_max)**2, N))**.5
    periods = np.exp(np.random.uniform(np.log(pmin), np.log(pmax), N))
    amps = np.exp(np.random.uniform(np.log(amp_min), np.log(amp_max), N))
    tau = np.exp(np.random.uniform(np.log(tau_min), np.log(tau_max),
                                   N))*periods
    return nspots, incl, periods, amps, tau


def generate_lcs(epic, N, nspot_min=50, nspot_max=500, incl_min=0,
                 incl_max=np.pi/4., amp_min=1e-8, amp_max=1e-4, pmin=.5,
                 pmax=90, tau_min=5, tau_max=20):
    """
    Generate N fake light curves based on the parameter limits.
    x (array): times of k2 lc.
    N (int): number of lcs.
    returns 2d array of light curves, the k2 light curve and a dictionary of
    true parameters.
    """

    x, y = k2lc(epic)

    params = [nspot_min, nspot_max, incl_min, incl_max, amp_min, amp_max,
              pmin, pmax, tau_min, tau_max]
    nspots, incl, periods, amps, tau = injection_params(N, params)
    true_params = {'nspots': nspots, 'incl': incl, 'periods': periods,
                   'amps': amps, 'tau': tau}

    xarr = np.zeros((len(x), N))
    yarr = np.zeros((len(x), N))
    for i in range(N):
        res0, res1 = mklc(x, nspot=nspots[i], incl=incl[i], tau=tau[i],
                          p=periods[i])
        med = np.median(res1[2, :])
        ys = (res1[2, :] / med - 1) * amps[i]  # median normalise and scale
        yarr[:, i] = ys + y
        xarr[:, i] = x
    if N == 1:
        return xarr.T[0], yarr.T[0], x, y, true_params
    return xarr, yarr, x, y, true_params

if __name__ == "__main__":
    xarr, yarr, x, y, true_params = generate_lcs("201131066", 1, amp_min=1,
                                                 amp_max=10)
    for i in range(len(xarr[0, :])):
        xs, ys = xarr[:, i], yarr[:, i]
        period, acf, lags, rvar = simple_acf(xs, ys)

        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(x, y, "b.")
        plt.plot(xs, ys, "k.", label="{0:.2f}".format(true_params["periods"][0]))
        plt.legend()
        plt.subplot(2, 1, 2)
        plt.plot(lags, acf)
        plt.axvline(period, color="r", label="{0:.2f}".format(period))
        plt.legend()
        plt.savefig("results/{0}".format(i))
