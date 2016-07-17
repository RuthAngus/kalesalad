from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import george
from george.kernels import ExpSine2Kernel, ExpSquaredKernel, WhiteKernel
import glob
import emcee
try:
    import corner
except:
    import triangle
import h5py
import subprocess
import scipy.optimize as spo
import time

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 20,
           'xtick.labelsize': 15,
           'ytick.labelsize': 15,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def lnprior(theta, plims):
    """
    plims is a tuple, list or array containing the lower and upper limits for
    the rotation period. These are logarithmic!
    theta = A, l, Gamma, s, P
    """
    if -20 < theta[0] < 20 and theta[4] < theta[1] and -20 < theta[2] < 20 \
    and -20 < theta[3] < 20 and plims[0] < theta[4] < plims[1]:
        return 0.
    return -np.inf


def Glnprior(theta, plims):
    """
    plims is a tuple, list or array containing the lower and upper limits for
    the rotation period. These are not logarithmic!
    """
    if -20 < theta[0] < 16 and -10 < theta[1] < 1 and -20 < theta[2] < 20 \
    and -20 < theta[3] < 20 and np.log(plims[0]) < theta[4] < np.log(plims[1]):
        return -.5 * ((theta[4] - np.log(plims[2]))/(.5*np.log(plims[2])))**2 \
                - .5 * ((theta[1] - .1*np.log(plims[2]))/(.1*plims[2]))**2
    return -np.inf

# lnprob
def lnprob(theta, x, y, yerr, plims):
    return lnlike(theta, x, y, yerr) + lnprior(theta, plims)

# lnlike
def lnlike(theta, x, y, yerr):
    theta = np.exp(theta)
    k = theta[0] * ExpSquaredKernel(theta[1]) \
            * ExpSine2Kernel(theta[2], theta[4]) + WhiteKernel(theta[3])
    gp = george.GP(k, solver=george.HODLRSolver)
    try:
        gp.compute(x, np.sqrt(theta[3]+yerr**2))
    except (ValueError, np.linalg.LinAlgError):
        return 10e25
    return gp.lnlikelihood(y, quiet=True)

# lnlike
def neglnlike(theta, x, y, yerr):
    theta = np.exp(theta)
    k = theta[0] * ExpSquaredKernel(theta[1]) \
            * ExpSine2Kernel(theta[2], theta[4])
    gp = george.GP(k)
    try:
        gp.compute(x, np.sqrt(theta[3]+yerr**2))
    except (ValueError, np.linalg.LinAlgError):
        return 10e25
    return -gp.lnlikelihood(y, quiet=True)

# make various plots
def make_plot(sampler, x, y, yerr, ID, DIR, traces=False, tri=False,
              prediction=True):

    nwalkers, nsteps, ndims = np.shape(sampler)
    flat = np.reshape(sampler, (nwalkers * nsteps, ndims))
    mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                      zip(*np.percentile(flat, [16, 50, 84], axis=0)))
    mcmc_result = np.array([i[0] for i in mcmc_result])
    print("\n", np.exp(np.array(mcmc_result[-1])), "period (days)", "\n")
    print(mcmc_result)
    np.savetxt("%s/%s_result.txt" % (DIR, ID), mcmc_result)

    fig_labels = ["$\ln(A)$", "$'ln(l)$", "$\ln(\Gamma)", "$\ln(\sigma)$",
                  "$P_{rot}$"]

    if traces:
        print("Plotting traces")
        for i in range(ndims):
            plt.clf()
            plt.plot(sampler[:, :, i].T, 'k-', alpha=0.3)
            plt.ylabel(fig_labels[i])
            plt.savefig("%s/%s_%s.png" % (DIR, ID, fig_labels[i]))

    if tri:
        print("Making triangle plot")
        flat[:, -1] = np.exp(flat[:, -1])
        try:
            fig = corner.corner(flat, labels=fig_labels)
        except:
            fig = triangle.corner(flat, labels=fig_labels)
        fig.savefig("%s/%s_triangle" % (DIR, ID))
        print("%s/%s_triangle.png" % (DIR, ID))

    if prediction:
        print("plotting prediction")
        theta = np.exp(np.array(mcmc_result))
        k = theta[0] * ExpSquaredKernel(theta[1]) \
                * ExpSine2Kernel(theta[2], theta[4])
        gp = george.GP(k, solver=george.HODLRSolver)
        gp.compute(x, yerr)
        xs = np.linspace(x[0], x[-1], 1000)
        mu, cov = gp.predict(y, xs)
        plt.clf()
        plt.errorbar(x-x[0], y, yerr=yerr, fmt="k.", capsize=0)
        plt.xlabel("$\mathrm{Time~(days)}$")
        plt.ylabel("$\mathrm{Normalised~Flux}$")
        plt.plot(xs, mu, color="CornFlowerBlue")
        plt.xlim(min(x), max(x))
        plt.savefig("%s/%s_prediction" % (DIR, ID))
        print("%s/%s_prediction.png" % (DIR, ID))

# take x, y, yerr and initial guess and do MCMC
def MCMC(theta_init, x, y, yerr, plims, burnin, run, ID, DIR, nwalkers=32,
         logsamp=True, plot_inits=False):

    # figure out whether x, y and yerr are arrays or lists of lists
    quarters = False
    if len(x) < 20:
        quarters = True

    print("\n", "log(theta_init) = ", theta_init)
    print("theta_init = ", np.exp(theta_init), "\n")

    if quarters:
        xl = [i for j in x for i in j]
        yl = [i for j in y for i in j]
        yerrl = [i for j in yerr for i in j]
        print("plotting inits")
        print(np.exp(theta_init))
        t = np.exp(theta_init)
        k = t[0] * ExpSquaredKernel(t[1]) * ExpSine2Kernel(t[2], t[3])
        gp = george.GP(k)
        gp.compute(xl, yerrl)
        xs = np.linspace(xl[0], xl[-1], 1000)
        mu, cov = gp.predict(yl, xs)

        plt.clf()
        plt.errorbar(xl, yl, yerr=yerrl, **reb)
        plt.plot(xs, mu, color="CornFlowerBlue")

        args = (xl, yl, yerrl)
        results = spo.minimize(neglnlike, theta_init, args=args)
        print("optimisation results = ", results.x)

        r = np.exp(results.x)
        k = r[0] * ExpSquaredKernel(r[1]) * ExpSine2Kernel(r[2], r[3])
        gp = george.GP(k)
        gp.compute(xl, yerrl)

        mu, cov = gp.predict(yl, xs)
        plt.plot(xs, mu, color="m", alpha=.5)
        plt.savefig("%s/%s_init" % (DIR, ID))
        print("%s/%s_init.png" % (DIR, ID))

    ndim, nwalkers = len(theta_init), nwalkers
    p0 = [theta_init+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
    args = (x, y, yerr, plims)

    lp = lnprob
    if quarters:  # if fitting each quarter separately, use a different lnprob
        lp = lnprob_split

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lp, args=args)

    print("burning in...")
    p0, lp, state = sampler.run_mcmc(p0, burnin)
    sampler.reset()
    print("production run...")
    p0, lp, state = sampler.run_mcmc(p0, run)

    # save samples
    f = h5py.File("%s/%s_samples.h5" % (DIR, ID), "w")
    data = f.create_dataset("samples", np.shape(sampler.chain))
    data[:, :] = np.array(sampler.chain)
    f.close()
    return sampler
