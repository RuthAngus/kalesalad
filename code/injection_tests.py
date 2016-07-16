import numpy as np
import matplotlib.pyplot as plt
from kalesalad import process_data
from mklc import mklc
from simple_acf import simple_acf
from Kepler_ACF import corr_run
import time


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


def generate_lcs(x, y, id, N, nspot_min=50, nspot_max=500, incl_min=0,
                 incl_max=np.pi/4., amp_min=1, amp_max=100, pmin=.5,
                 pmax=90, tau_min=5, tau_max=20):
    """
    Generate N fake light curves based on the parameter limits.
    x (array): times of k2 lc.
    N (int): number of lcs.
    returns 2d array of light curves, the k2 light curve and a dictionary of
    true parameters.
    """

    params = [nspot_min, nspot_max, incl_min, incl_max, amp_min*rvar,
              amp_max*rvar, pmin, pmax, tau_min, tau_max]
    nspots, incl, periods, amps, tau = injection_params(N, params)
    true_params = {'nspots': nspots, 'incl': incl, 'periods': periods,
                   'amps': amps, 'tau': tau}

    xarr = np.zeros((len(x), N))
    yarr = np.zeros((len(x), N))
    for i in range(N):
        print(i, "of", N)
        res0, res1 = mklc(x, nspot=nspots[i], incl=incl[i], tau=tau[i],
                          p=periods[i])
        med = np.median(res1[2, :])
        ys = (res1[2, :] / med - 1) * amps[i]  # median normalise and scale
        yarr[:, i] = ys + y
        xarr[:, i] = x
    if N == 1:
        return xarr.T[0], yarr.T[0], x, y, true_params

    # save the results
    np.savetxt("lcs.txt", yarr.T)
    np.savetxt("xs.txt", xarr.T)
    np.savetxt("truth.txt", np.vstack((nspots, incl, periods, amps, tau)).T)

    return xarr, yarr, true_params

if __name__ == "__main__":
    N = 100  # number of light curves to simulate
    nspot_min = 50  # minimum number of spots
    nspot_max = 500  # maximum number of spots
    incl_min = 0  # minimum inclination
    incl_max = np.pi/4.  # maximum inclination
    amp_min = 1  # minimum amplitude (multiple of range of variability)
    amp_max = 100  # maximum amplitude (see above)
    pmin = .5  # minimum period (days)
    pmax = 90  # maximum period (days)
    tau_min = 5  # minimum spot lifetime (multiple of rotation period)
    tau_max = 20  # maximum spot lifetime (see above)

    # load k2 lc
    epic = "201131066"
    x, y = k2lc(epic)
    period, acf, lags, rvar, peaks, dips, leftdips, rightdips, bigpeaks \
        = simple_acf(x, y)

    start = time.time()
    acf, lags, period, err, locheight = corr_run(x, y)
    end = time.time()
    print("time = ", end-start)

    # xarr, yarr, true_params = generate_lcs(x, y, epic, N,
    #                                              nspot_min=nspot_min,
    #                                              nspot_max=nspot_max,
    #                                              incl_min=incl_min,
    #                                              incl_max=incl_max,
    #                                              amp_min=amp_min,
    #                                              amp_max=amp_max,
    #                                              pmin=pmin, pmax=pmax,
    #                                              tau_min=tau_min,
    #                                              tau_max=tau_max)
    # periods = true_params["periods"]
    # amps = true_params["amps"]

    # load simulations
    yarr = np.genfromtxt("lcs.txt").T
    xarr = np.genfromtxt("xs.txt").T
    _, _, periods, amps, _ = np.genfromtxt("truth.txt").T

    # recover and make plots
    recovered = []
    for i in range(N):
        xs, ys = xarr[:, i], yarr[:, i]
        acf, lags, period, err, locheight = corr_run(xs, ys)
        recovered.append(period)

        print(i, "of", N)
        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(x, y, "b.")
        plt.plot(xs, ys, "k.",
                 label="{0:.2f}".format(periods[i]))
        plt.legend()
        plt.subplot(2, 1, 2)
        plt.plot(lags, acf)
        plt.axvline(period, color="r", label="{0:.2f}".format(locheight))
        plt.legend()
        plt.savefig("results/simulations/{0}".format(i))

    plt.clf()
    periods, recovered = np.array(periods[:N]), np.array(recovered)
    m = periods > 0
    plt.scatter(periods[m], recovered[m], marker="o", c=amps[m],
                edgecolor="")
    plt.colorbar()
    plt.ylim(0, 100)
    xs = np.linspace(min(periods), max(periods), 100)
    plt.plot(xs, xs, "k--")
    plt.plot(xs, .5 * xs, "k--")
    plt.savefig("test")

    success = periods[np.abs(periods-recovered) < .1*periods]
    print(len(success) / float(N), "success fraction")
