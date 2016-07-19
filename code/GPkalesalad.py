from __future__ import print_function
import numpy as np
from GProtation import MCMC, make_plot, lnprob
from kalesalad import process_data
import os
from simple_acf import simple_acf
import time
import emcee
import h5py
import matplotlib.pyplot as plt


def recover_injections(id, x, y, yerr, fn, burnin, run, interval, tol,
                       npts=10, nwalkers=32, plot=True):
    """
    Take x, y, yerr, calculate ACF period for initialisation and do MCMC.
    npts: number of points per period.
    """

    p_init, acf_smooth, lags, _, _, _, _, _, _ = simple_acf(x, y)

    print("acf period = ", p_init)

    if p_init < .1:  # prevent unphysical periods
            p_init = 10.

    # Format data
    plims = np.log([p_init - tol*p_init, p_init + tol*p_init])

    print(p_init, np.exp(plims))

    sub = int(p_init / float(npts) * 48)  # 10 points per period
    ppd = 48. / sub
    ppp = ppd * p_init
    print("sub = ", sub, "points per day =", ppd, "points per period =", ppp)
    # subsample
    xsub, ysub, yerrsub = x[::sub], y[::sub], yerr[::sub]
    xb, yb, yerrb = x, y, yerr
    xb, yb, yerrb = x[:100], y[:100], yerr[:100]
    plt.clf()
    plt.plot(xb, yb, "k.")
    plt.savefig("gptest")

    theta_init = np.log([np.exp(-5), np.exp(7), np.exp(.6), np.exp(-16),
                         p_init])

    print("\n", "log(theta_init) = ", theta_init)
    print("theta_init = ", np.exp(theta_init), "\n")

    # set up MCMC
    ndim, nwalkers = len(theta_init), nwalkers
    p0 = [theta_init+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
    args = (xb, yb, yerrb, plims)
    lp = lnprob

    # time the lhf call
    start = time.time()
    print("lnprob = ", lp(theta_init, xb, yb, yerrb, plims))
    end = time.time()
    tm = end - start
    print("1 lhf call takes ", tm, "seconds")
    print("burn in will take", tm * nwalkers * burnin, "s")
    print("run will take", tm * nwalkers * run, "s")
    print("total = ", (tm * nwalkers * run + tm * nwalkers * burnin)/60,
          "mins")

    # run MCMC
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lp, args=args)
    print("burning in...")
    start = time.time()
    p0, lp, state = sampler.run_mcmc(p0, burnin)
    sampler.reset()
    print("production run...")
    p0, lp, state = sampler.run_mcmc(p0, run)
    end = time.time()
    print("actual time = ", end - start)

    # save samples
    f = h5py.File("%s_samples.h5" % id, "w")
    data = f.create_dataset("samples", np.shape(sampler.chain))
    data[:, :] = np.array(sampler.chain)
    f.close()

    # make various plots
    if plot:
        with h5py.File("%s_samples.h5" % id, "r") as f:
            samples = f["samples"][...]
        mcmc_result = make_plot(samples, xsub, ysub, yerrsub, id, fn,
                                traces=True, tri=True, prediction=True)

if __name__ == "__main__":

    c = "01"
    epic = "201131793"
    path = "data/c01/201100000/31793"
    end = "kepler_v1.0_lc.fits"
    file = "{0}/hlsp_everest_k2_llc_{1}-c{2}_{3}".format(path, epic, c, end)
    if os.path.exists(file):
        x, y = process_data(file)

        burnin, run, npts, tol = 1000, 100000, 10, .4  # MCMC. max npts is 48
        yerr = np.ones_like(y) * 1e-5
        interval = 0.02043365  # assume for long cadence
        recover_injections(id, x, y, yerr, path, burnin, run, interval, tol,
                           npts, nwalkers=12, plot=True)
    else:
        print(file, "file not found")
