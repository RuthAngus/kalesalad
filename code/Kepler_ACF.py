# This version of the ACF code is designed to run without the stupid index step
# It is also totally stripped down and only saves the information that I actually use.
# the corr_run function takes

import scipy
from numpy.random import normal
import matplotlib.image as mpimg
import random
import numpy as np
import atpy
import pylab
import copy
import glob
# from sets import Set
import collections
no_rpy = True
import scipy.io
from scipy import signal
import KOI_tools_b12 as kt
import filter
import gls
import mpfit
import pyfits
import matplotlib.pyplot as pl

gap_days = 0.02043365  # assume for long cadence
jump_arr = scipy.array([131.51139, 169.51883, 169.75000, 182.00000, 200.31000,
                       231.00000, 246.19000, 256.00000, 260.22354, 281.00000,
                       291.00000, 322.00000, 352.37648, 373.23000, 384.00000,
                       398.00000, 443.48992, 475.50000, 504.00000, 539.44868,
                       567.00000, 599.00000, 630.17387, 661.00000, 691.00000,
                       711.20000, 735.36319, 762.00000, 808.51558, 845.00000,
                       874.50000, 906.84469, 937.00000, 970.00000, 1001.20718,
                       1032.50000, 1063.50000 ,1071.00000, 1093.60000])

def corr_run(time, flux):

    id_list = [0]
    # Create empty arrays
    acf_peak_per = scipy.ones(len(id_list)) * -9999.0
    sine_per = scipy.ones(len(id_list)) * -9999.0
    sine_height = scipy.ones(len(id_list)) * -9999.0
    med_dlag_per = scipy.ones(len(id_list)) * -9999.0
    dlag_per_err = scipy.ones(len(id_list)) * -9999.0
    h1 = scipy.ones(len(id_list)) * -9999.0
    hlocgrad = scipy.ones(len(id_list)) * -9999.0
    hloc_grad_scatter = scipy.ones(len(id_list)) * -9999.0
    width_grad = scipy.ones(len(id_list)) * -9999.0
    width_grad_scatter = scipy.ones(len(id_list)) * -9999.0
    w1 = scipy.ones(len(id_list)) * -9999.0
    lh1 = scipy.ones(len(id_list)) * -9999.0
    num_of_peaks = scipy.ones(len(id_list)) * -9999.0
    harmonic_det = scipy.ones(len(id_list)) * -9999.0
    amp_all = scipy.ones(len(id_list)) * -9999.0
    amp_per = scipy.ones(len(id_list)) * -9999.0
    period = scipy.ones(len(id_list)) * -9999.0
    mdn = np.median(flux)
    flux = flux-mdn
    lc_tab = atpy.Table()
    lc_tab.add_column('time', time)
    lc_tab.add_column('flux', flux)
    lc_tab.add_column('flux_pdc', flux)

    qt_max = [0., max(lc_tab.time)]
    tablen = 1
    x = 0

    # max period searched for is len(flux) / 2
    max_psearch_len = len(lc_tab.flux) / 2.0

    # Calculate ACF
    print('Calculating ACF...')

    acf_tab, acf_per_pos, acf_per_height, acf_per_err, locheight, asym,  = \
        acf_calc(time = lc_tab.time, flux = lc_tab.flux, interval = gap_days, \
                 kid = x, max_psearch_len = max_psearch_len)
    acf_smooth, lags = acf_tab["acf_smooth"], acf_tab["lags_days"]

    pgram_tab, sine_per[x], sine_height[x] = \
        pgram_calc(time = lc_tab.time, flux = lc_tab.flux, \
                   interval = gap_days, kid = x, max_psearch_len = max_psearch_len)

    # plot and calculate acf peak statistics
    if acf_per_pos[0] != -9999:
        med_dlag_per[x], dlag_per_err[x], acf_peak_per[x], h1[x], w1[x], lh1[x], \
            hlocgrad[x], hloc_grad_scatter[x], width_grad[x], width_grad_scatter[x], \
            num_of_peaks[x], harmonic_det[x], sel_peaks, one_peak_only, peak_ratio =\
            plot_stats(lc_tab.time, lc_tab.flux, x, acf_per_pos, \
                       acf_per_height, acf_per_err, locheight, asym)

        n_s = 'k'

        period[x] = acf_peak_per[x]
        locheight = locheight[x]
        if locheight == [-9999]:
            locheight = 0
        if locheight < 0.1:
            period[x] = 0

        # variability stats
        amp_all[x], amp_per[x], per_cent, var_arr_real = \
            calc_var(kid = x, time_in = lc_tab.time, \
                     flux = lc_tab.flux, period = acf_peak_per[x])

        maxpts = 40.0
        if scipy.floor(lc_tab.time.max() / acf_peak_per[x]) < maxpts:
            maxpts = float(scipy.floor(lc_tab.time.max() / acf_peak_per[x]))
        inc = lc_tab.time - lc_tab.time.min() <= (maxpts*acf_peak_per[x])

        t_stats = atpy.Table()
        t_stats.add_column('acf_per_pos', acf_per_pos)
        t_stats.add_column('acf_per_height', acf_per_height)
        t_stats.add_column('acf_per_err', acf_per_err)
        t_stats.add_column('asym', asym)
        t_stats.add_column('locheight', locheight)

        if dlag_per_err[x] == 0.:
             error = acf_per_err[x]
        else: error = dlag_per_err[x]

    else:
        blank = np.array([0,0])

    t = atpy.Table()
    t.add_column('period', period) #period
    t.add_column('sine_per', sine_per) #sine period
    t.add_column('sine_height', sine_height)
    t.add_column('acf_peak_per', acf_peak_per)
    t.add_column('med_dlag_per', med_dlag_per)
    t.add_column('dlag_per_err', dlag_per_err) #error
    t.add_column('h1', h1)
    t.add_column('w1', w1)
    t.add_column('lh1', lh1)
    t.add_column('hlocgrad', hlocgrad)
    t.add_column('hloc_grad_scatter', hloc_grad_scatter)
    t.add_column('width_grad', width_grad)
    t.add_column('width_grad_scatter', width_grad_scatter)
    t.add_column('num_of_peaks', num_of_peaks)
    t.add_column('harmonic_det', harmonic_det)
    t.add_column('amp_all', amp_all)
    t.add_column('amp_per', amp_per)
    return acf_smooth, lags, period[0], dlag_per_err, locheight

def acf_calc(time, flux, interval, kid, max_psearch_len):

    ''' Calculate ACF, calls error calc function'''
    lags, acf, lines, axis = pylab.acorr(flux, maxlags = max_psearch_len)

    #convolve smoothing window with Gaussian kernel
    gauss_func = lambda x,sig: 1./np.sqrt(2*np.pi*sig**2) * \
                 np.exp(-0.5*(x**2)/(sig**2)) #define a Gaussian
    #create the smoothing kernel
    conv_func = gauss_func(np.arange(-28,28,1.),9.)

    acf_smooth = np.convolve(acf,conv_func,mode='same') #and convolve
    lenlag = len(lags)
    lags = lags[int(lenlag/2.0):lenlag][:-1] * interval
    acf = acf[int(lenlag/2.0): lenlag][0:-1]
    acf_smooth = acf_smooth[int(lenlag/2.0): lenlag][1:]

    # find max using usmoothed acf (for plot only)
    max_ind_us, max_val_us = extrema(acf, max = True, min = False)

    # find max/min using smoothed acf
    max_ind_s, max_val_s = extrema(acf_smooth, max = True, min = False)
    min_ind_s, min_val_s = extrema(acf_smooth, max = False, min = True)
    maxmin_ind_s, maxmin_val_s = extrema(acf_smooth, max = True, min = True)

    if len(max_ind_s) > 0 and len(min_ind_s) > 0:
        # ensure no duplicate peaks are detected
        t_max_s = atpy.Table()
        t_max_s.add_column('ind', max_ind_s)
        t_max_s.add_column('val', max_val_s)
        t_min_s = atpy.Table()
        t_min_s.add_column('ind', min_ind_s)
        t_min_s.add_column('val', min_val_s)
        t_maxmin_s = atpy.Table()
        t_maxmin_s.add_column('ind', maxmin_ind_s)
        t_maxmin_s.add_column('val', maxmin_val_s)

        ma_i = collections.Counter(t_max_s.ind)
        dup_arr = [i for i in ma_i if ma_i[i]>1]
        if len(dup_arr) > 0:
            for j in scipy.arange(len(dup_arr)):
                tin = t_max_s.where(t_max_s.ind != dup_arr[j])
                tout = t_max_s.where(t_max_s.ind == dup_arr[j])
                tout = tout.rows([0])
                tin.append(tout)
            t_max_s = copy.deepcopy(tin)

        ma_i = collections.Counter(t_min_s.ind)
        dup_arr = [i for i in ma_i if ma_i[i]>1]
        if len(dup_arr) > 0:
            for j in scipy.arange(len(dup_arr)):
                tin = t_min_s.where(t_min_s.ind != dup_arr[j])
                tout = t_min_s.where(t_min_s.ind == dup_arr[j])
                tout = tout.rows([0])
                tin.append(tout)
            t_min_s = copy.deepcopy(tin)

        ma_i = collections.Counter(t_maxmin_s.ind)
        dup_arr = [i for i in ma_i if ma_i[i]>1]
        if len(dup_arr) > 0:
            for j in scipy.arange(len(dup_arr)):
                tin = t_maxmin_s.where(t_maxmin_s.ind != dup_arr[j])
                tout = t_maxmin_s.where(t_maxmin_s.ind == dup_arr[j])
                tout = tout.rows([0])
                tin.append(tout)
            t_maxmin_s = copy.deepcopy(tin)

        t_max_s.sort('ind')
        t_min_s.sort('ind')
        t_maxmin_s.sort('ind')

        # relate max inds to lags
        maxnum = len(t_max_s.ind)
        acf_per_pos = lags[t_max_s.ind]
        acf_per_height = acf[t_max_s.ind]

        # Calculate peak widths, asymmetries etc
        acf_per_err, locheight, asym= \
            calc_err(kid = kid, lags = lags, acf = acf, inds = \
            t_maxmin_s.ind, vals = t_maxmin_s.val, maxnum = maxnum)

    else:
        acf_per_pos = scipy.array([-9999])
        acf_per_height = scipy.array([-9999])
        acf_per_err = scipy.array([-9999])
        locheight = scipy.array([-9999])
        asym = scipy.array([-9999])

    # save corrected LC and ACF
    t_lc = atpy.Table()
    t_lc.add_column('time', time)
    t_lc.add_column('flux', flux)

    t_acf = atpy.Table()
    t_acf.add_column('lags_days', lags)
    t_acf.add_column('acf', acf)
    t_acf.add_column('acf_smooth', acf_smooth)

    return t_acf, acf_per_pos, acf_per_height, acf_per_err, locheight, asym

def pgram_calc(time, flux, interval, kid, max_psearch_len):
    ''' Calculate Sine Fitting Periodogram'''
    # Calculate sine lsq periodogram
    pmin = 0.1
    pmax = interval * max_psearch_len
    nf = 1000
    sinout = gls.sinefit(time, flux, err = None, pmin = pmin, pmax = pmax, nper = nf, \
                             doplot = False, return_periodogram = True)

    sine_per = sinout[0]
    sine_height = max(sinout[5])

    # save periodogram
    t_pg = atpy.Table()
    t_pg.add_column('period', sinout[4])
    t_pg.add_column('pgram', sinout[5])

    return t_pg, sine_per, sine_height

def calc_err(kid, lags, acf, inds, vals, maxnum):
    ''' Calculate peak widths, heights and asymmetries '''
    if len(inds) == 0: return -9999, -9999, -9999
    acf_per_err = scipy.ones(maxnum) * -9999
    asym = scipy.ones(maxnum) * -9999
    mean_height = scipy.ones(maxnum) * -9999

    acf_ind = scipy.r_[0:len(acf):1]
    num = len(vals)
    if maxnum*2 > num: maxnum -= 1
    # loop through maxima, assuming 1st index is for minima
    for i in scipy.arange(maxnum):
        # find values and indices of centre left and right
        centre_v = vals[2*i+1]
        centre_i = inds[2*i+1]

        # select value half way between max and min to calc width and asymmetry
        if 2*i + 2 >= num:
            # if it goes beyond limit of lags
            acf_per_err[i] = -9999
            mean_height[i] = -9999
            asym[i] = -9999
        else:
            left_v = vals[2*i]
            left_i = inds[2*i]
            right_v = vals[2*i+2]
            right_i = inds[2*i+2]

            sect_left_acf = acf[left_i:centre_i+1]
            sect_left_ind = acf_ind[left_i:centre_i+1]
            sect_right_acf = acf[centre_i:right_i+1]
            sect_right_ind = acf_ind[centre_i:right_i+1]

            height_r = centre_v - right_v
            height_l = centre_v - left_v

            # calc height from peak down 0.5 * mean of side heights
            mean_height[i] = 0.5*(height_l + height_r)
            mid_height = centre_v - 0.5*mean_height[i]
            if mid_height <= min(sect_left_acf): mid_height = min(sect_left_acf)
            if mid_height <= min(sect_right_acf): mid_height = min(sect_right_acf)

            sect_left_acf_r = sect_left_acf[::-1]
            sect_left_ind_r = sect_left_ind[::-1]
            for j in scipy.arange(len(sect_left_acf)):
                if sect_left_acf_r[j] <= mid_height:
                    if j == 0:
                        lag_mid_left = lags[sect_left_ind_r[j]]
                        break
                    else:
                        pt1 = sect_left_acf_r[j-1]
                        pt2 = sect_left_acf_r[j]
                        lag1 = lags[sect_left_ind_r[j-1]]
                        lag2 = lags[sect_left_ind_r[j]]
                        if pt1 < pt2:
                            ptarr = scipy.array([pt1, pt2])
                            lagarr = scipy.array([lag1, lag2])
                        else:
                            ptarr = scipy.array([pt2, pt1])
                            lagarr = scipy.array([lag2, lag1])
                        f_left = scipy.interpolate.interp1d(ptarr, lagarr)
                        lag_mid_left = f_left(mid_height)
                        break

            for j in scipy.arange(len(sect_right_acf)):
                if sect_right_acf[j] <= mid_height:
                    if j == 0:
                        lag_mid_right = lags[sect_right_ind[j]]
                        break
                    else:
                        pt1 = sect_right_acf[j-1]
                        pt2 = sect_right_acf[j]
                        lag1 = lags[sect_right_ind[j-1]]
                        lag2 = lags[sect_right_ind[j]]
                        if pt1 < pt2:
                            ptarr = scipy.array([pt1, pt2])
                            lagarr = scipy.array([lag1, lag2])
                        else:
                            ptarr = scipy.array([pt2, pt1])
                            lagarr = scipy.array([lag2, lag1])
                        f_right = scipy.interpolate.interp1d(ptarr, lagarr)
                        lag_mid_right = f_right(mid_height)
                        break

            pos_l = lag_mid_right - lags[centre_i]
            pos_r = lags[centre_i] - lag_mid_left
            asym[i] = pos_r / pos_l
            acf_per_err[i] = pos_l + pos_r
            if asym[i] <= 0:
                acf_per_err[i] = -9999
                mean_height[i] = -9999
                asym[i] = -9999

    return acf_per_err, mean_height, asym


def plot_stats(time, flux, kid_x, acf_per_pos_in, acf_per_height_in, acf_per_err_in, locheight_in, asym_in):
    ''' Plot and calculate statistics of peaks '''

    acf_per_pos_in = acf_per_pos_in[asym_in != -9999]
    acf_per_height_in = acf_per_height_in[asym_in != -9999]
    acf_per_err_in = acf_per_err_in[asym_in != -9999]
    locheight_in = locheight_in[asym_in != -9999]
    asym_in = asym_in[asym_in != -9999]

    if len(acf_per_pos_in) == 0: return -9999, -9999, -9999, -9999, -9999, -9999,\
        -9999, -9999, -9999, -9999, 0, -9999, -9999, 1, 0.0

    x = 10  #number of periods used for calc
    hdet = 0  #start with 0 harmonic, set to 1 if 1/2P is 1st peak
    # deal with cases where 1/2 P is 1st peak
    if len(acf_per_pos_in) >= 2:
        one_peak_only = 0
        ind = scipy.r_[1:len(acf_per_pos_in)+1:1]
        if locheight_in[1] > locheight_in[0]:
            hdet = 1  # mark harmonic found
            pk1 = acf_per_pos_in[1]
            acf_per_pos_in = acf_per_pos_in[1:]
            acf_per_height_in = acf_per_height_in[1:]
            acf_per_err_in = acf_per_err_in[1:]
            locheight_in = locheight_in[1:]
            asym_in = asym_in[1:]
            '''if 1 == 1:
            pknumin = int(raw_input('pk in: '))
            if pknumin > 0: hdet = 1  # mark harmonic found
            pk1 = acf_per_pos_in[pknumin]
            acf_per_pos_in = acf_per_pos_in[pknumin:]
            acf_per_height_in = acf_per_height_in[pknumin:]
            acf_per_err_in = acf_per_err_in[pknumin:]
            locheight_in = locheight_in[pknumin:]
            asym_in = asym_in[pknumin:]'''
        else:
            pk1 = acf_per_pos_in[0]
    else:
        one_peak_only = 1
        pk1 = acf_per_pos_in[0]


    # select only peaks which are ~multiples of 1st peak (within phase 0.2)
    acf_per_pos_0_test = scipy.append(0,acf_per_pos_in)
    ind_keep = scipy.r_[0:len(acf_per_pos_0_test):1]
    fin = False
    while fin == False:
        delta_lag_test = acf_per_pos_0_test[1:] - acf_per_pos_0_test[:-1]
        delta_lag_test = scipy.append(0, delta_lag_test)
        phase = ((delta_lag_test % pk1) / pk1)
        phase[phase > 0.5] -= 1.0
        excl = abs(phase) > 0.2
        ind_temp = scipy.r_[0:len(delta_lag_test):1]
        if len(phase[excl]) == 0: break
        else:
            ind_rem = ind_temp[excl][0]
            rem = acf_per_pos_0_test[ind_rem]
            ind_keep = ind_keep[acf_per_pos_0_test != rem]
            acf_per_pos_0_test = acf_per_pos_0_test[acf_per_pos_0_test != rem]

    ind_keep = ind_keep[1:] - 1
    keep_pos = acf_per_pos_in[ind_keep]
    keep_pos = scipy.append(0, keep_pos)
    delta_keep_pos = keep_pos[1:] - keep_pos[:-1]
    # remove very small delta lags points (de-noise peak detections)
    if len(ind_keep) > 1:
        ind_keep = ind_keep[delta_keep_pos > 0.3*delta_keep_pos[0]]
        delta_keep_pos = delta_keep_pos[delta_keep_pos > 0.3*delta_keep_pos[0]]

    if len(ind_keep) > 1:
        ind_gap = ind_keep[delta_keep_pos > 2.2*delta_keep_pos[0]]
        if len(ind_gap) != 0:
            ind_keep = ind_keep[ind_keep < ind_gap[0]]
            delta_keep_pos = delta_keep_pos[ind_keep < ind_gap[0]]

    # limit to x lags for plot and calc
    if len(acf_per_pos_in[ind_keep]) < x: x = len(acf_per_pos_in[ind_keep])

    acf_per_pos = acf_per_pos_in[ind_keep][:x]
    acf_per_height = acf_per_height_in[ind_keep][:x]
    acf_per_err = acf_per_err_in[ind_keep][:x]
    asym = asym_in[ind_keep][:x]
    locheight = locheight_in[ind_keep][:x]

    if len(acf_per_pos) == 1:
        return -9999, 0.0, pk1, acf_per_height[0], acf_per_err[0], locheight[0], -9999, -9999, -9999, -9999, 1, hdet, -9999, 1, 0.0

    ''' Delta Lag '''
    acf_per_pos_0 = scipy.append(0,acf_per_pos)
    delta_lag = acf_per_pos_0[1:] - acf_per_pos_0[:-1]
    av_delt = scipy.median(delta_lag)
    delt_mad = 1.483*scipy.median(abs(delta_lag - av_delt)) # calc MAD
    delt_mad = delt_mad / scipy.sqrt(float(len(delta_lag)-1.0)) # err = MAD / sqrt(n-1)
    med_per = av_delt
    mad_per_err = delt_mad

    if len(acf_per_pos) > 2 and no_rpy == False:
        st_line = lambda p, x: p[0] + p[1] * x
        '''st_line_err = lambda p, x, y, fjac: [0, (y - st_line(p, x)), None]
        fa = {'x': acf_per_pos, 'y': locheight}
        p = [scipy.median(locheight), 0.0]
        m = mpfit.mpfit(st_line_err, p, functkw = fa, quiet = True)
        p = m.params
        errs = m.perror
        h_grad = p[1]
        h_timescale = 1.0 / h_grad
        h_grad_err = errs[1]'''
        r.assign('xdf1', acf_per_pos)
        r.assign('ydf1', locheight)
        p1 = r('''
        xdf <- c(xdf1)
        ydf <- c(ydf1)
        library(quantreg)
        rqmodel1 <- rq(ydf~xdf)
        plot(xdf,ydf)
        abline(rqmodel1,col=3)
        sumr = summary(rqmodel1, se = 'boot')
        grad <- rqmodel1[['coefficients']][['xdf']]
        interc <- rqmodel1[['coefficients']][['(Intercept)']]
        err_grad <- sumr[[3]][[3]]
        err_interc <- sumr[[3]][[4]]
        resids = resid(rqmodel1)
        output <- c(resids, grad, interc, err_grad, err_interc)
        ''')
        res =scipy.array(p1[0:len(acf_per_pos)])
        pqr = scipy.array(p1[len(acf_per_pos):])
        h_grad = pqr[0]
        h_timescale = 1.0 / h_grad
        h_grad_err = pqr[3]
        h_grad_scatter = sum((st_line([pqr[1],pqr[0]], acf_per_pos) - locheight) ** 2) / scipy.sqrt(float(len(acf_per_pos)))
    else:
        h_grad = -9999
        h_timescale = -9999
        h_grad_err = -9999
        h_grad_scatter = -9999

    if len(acf_per_pos) > 2 and no_rpy == False:
        r.assign('xdf1', acf_per_pos)
        r.assign('ydf1', acf_per_err)
        p1 = r('''
        xdf <- c(xdf1)
        ydf <- c(ydf1)
        library(quantreg)
        rqmodel1 <- rq(ydf~xdf)
        plot(xdf,ydf)
        abline(rqmodel1,col=3)
        sumr = summary(rqmodel1, se = 'boot')
        grad <- rqmodel1[['coefficients']][['xdf']]
        interc <- rqmodel1[['coefficients']][['(Intercept)']]
        err_grad <- sumr[[3]][[3]]
        err_interc <- sumr[[3]][[4]]
        resids = resid(rqmodel1)
        output <- c(resids, grad, interc, err_grad, err_interc)
        ''')
        res =scipy.array(p1[0:len(acf_per_pos)])
        pqr = scipy.array(p1[len(acf_per_pos):])
        w_grad = pqr[0]
        w_timescale = 1.0 / w_grad
        w_grad_err = pqr[3]
        w_grad_scatter = sum((st_line([pqr[1],pqr[0]], acf_per_pos) - acf_per_err) ** 2) / scipy.sqrt(float(len(acf_per_pos)))
    else:
        w_grad = -9999
        w_timescale = -9999
        w_grad_err = -9999
        w_grad_scatter = -9999

    peak_ratio = len(acf_per_pos)/len(acf_per_pos_in)

    return med_per, mad_per_err, pk1, acf_per_height[acf_per_pos == pk1][0], \
            acf_per_err[acf_per_pos == pk1][0], \
            locheight[acf_per_pos == pk1][0], h_grad, h_grad_scatter, w_grad, \
            w_grad_scatter, len(acf_per_pos), hdet, acf_per_pos, \
            one_peak_only, peak_ratio

def calc_var(kid = None, time_in = None, flux = None, period = None):
    ''' calculate 5-95th percentile of flux (i.e. amplitude) whole LC and in
    each period block '''

    # Of whole LC...
    sort_flc= sorted(flux)
    fi_ind_flc = int(len(sort_flc) * 0.05)
    nifi_ind_flc = int(len(sort_flc) * 0.95)
    amp_flc = sort_flc[nifi_ind_flc] - sort_flc[fi_ind_flc]

    # Of each block...
    if period > 0:
        num = int(scipy.floor(len(time_in) / period))
        var_arr = scipy.zeros(num) + scipy.nan
        per_cent = scipy.zeros(num) + scipy.nan
        for i in scipy.arange(num-1):
            t_block = time_in[(time_in-time_in.min() >= i*period) * \
                    (time_in-time_in.min() < (i+1)*period)]
            f_block = flux[(time_in-time_in.min() >= i*period) * \
                    (time_in-time_in.min() < (i+1)*period)]
            if len(t_block) < 2: continue
            sort_f_block = sorted(f_block)
            fi_ind = int(len(sort_f_block) * 0.05)
            nifi_ind = int(len(sort_f_block) * 0.95)
            range_f = sort_f_block[nifi_ind] - sort_f_block[fi_ind]
            var_arr[i] = range_f
            per_cent[i] = scipy.mean(t_block)

        var_arr_real = var_arr[scipy.isfinite(var_arr) == True]
        per_cent = per_cent[scipy.isfinite(var_arr) == True]
        var_med = scipy.median(var_arr_real)

    else:
        var_med = -9999
        per_cent = scipy.array([-9999])
        var_arr_real = scipy.array([-9999])

    return amp_flc, var_med, per_cent, var_arr_real

def extrema(x, max = True, min = True, strict = False, withend = False):
	"""
	This function will index the extrema of a given array x.
	Options:
		max		If true, will index maxima
		min		If true, will index minima
		strict		If true, will not index changes to zero gradient
		withend	If true, always include x[0] and x[-1]
	This function will return a tuple of extrema indexies and values
	"""
	# This is the gradient
	from numpy import zeros
	dx = zeros(len(x))
	from numpy import diff
	dx[1:] = diff(x)
	dx[0] = dx[1]

	# Clean up the gradient in order to pick out any change of sign
	from numpy import sign
	dx = sign(dx)

	# define the threshold for whether to pick out changes to zero gradient
	threshold = 0
	if strict:
		threshold = 1

	# Second order diff to pick out the spikes
	d2x = diff(dx)

	if max and min:
		d2x = abs(d2x)
	elif max:
		d2x = -d2x

	# Take care of the two ends
	if withend:
		d2x[0] = 2
		d2x[-1] = 2

	# Sift out the list of extremas
	from numpy import nonzero
	ind = nonzero(d2x > threshold)[0]
	return ind, x[ind]
