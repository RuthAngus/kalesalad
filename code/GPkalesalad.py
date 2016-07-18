import numpy as np
from GProtation import MCMC, make_plot, lnprob
from kalesalad import process_data
import os
from Kepler_ACF import corr_run
import time
import emcee
import h5py


def recover_injections(id, x, y, yerr, fn, burnin, run, interval, tol,
                       npts=10, nwalkers=32, plot=True):
    """
    Take x, y, yerr, calculate ACF period for initialisation and do MCMC.
    npts: number of points per period.
    """

    acf, lags, p_init, err, locheight = corr_run(x, y)
    print("acf period, err = ", p_init)

    if p_init < .1:  # prevent unphysical periods
            p_init = 1.

    # Format data
    plims = np.log([p_init - tol*p_init, p_init + tol*p_init])

    print(p_init, np.exp(plims))

    sub = int(p_init / float(npts) * 48)  # 10 points per period
    print(sub, p_init, npts)
    assert 0
    ppd = 48. / sub
    ppp = ppd * p_init
    print("sub = ", sub, "points per day =", ppd, "points per period =", ppp)
    # subsample
    xsub, ysub, yerrsub = x[::sub], y[::sub], yerr[::sub]
    xb, yb, yerrb = x, y, yerr

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

        burnin, run, npts, tol = 500, 1000, 50, .4  # MCMC. max npts is 48
        yerr = np.ones_like(y) * 1e-5
        interval = 0.02043365  # assume for long cadence
        recover_injections(id, x, y, yerr, path, burnin, run, interval, tol,
                           npts, nwalkers=12, plot=True)
        assert 0
