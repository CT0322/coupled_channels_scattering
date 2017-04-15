################################################################################
#
# Author: Bastian Knippschild (b.knippschild@gmx.de)
# Date:   Februar 2015
#
# Copyright (C) 2015 Bastian Knippschild
# 
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with tmLQCD. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
#
# Function: Functions to fit and plot.
#
# For informations on input parameters see the description of the function.
#
################################################################################

import sys
from scipy.optimize import leastsq
import scipy.stats
import numpy as np
import analyze_fcts as af
from .plot import *

def constrain(p):
	a, b = 0., 0.
	if p[0] > 80. or p[0] < 0.:
		a = 10000.
	if p[1] > 1.0 or p[1] < 0.01:
		b = 10000.
	return a+b

def fitting(fitfunc, X, Y, start_parm, correlated=True, verbose=False):
    """A function that fits a correlation function.

    This function fits the given function fitfunc to the data given in X and Y.
    The function needs some start values, given in start_parm, and can use a
    correlated or an uncorrelated fit.

    Args:
        fitfunc: The function to fit to the data.
        X: The time slices.
        Y: The bootstrap samples of the data.
        start_parm: The starting parameters for the fit.
        correlated: Flag to use a correlated or uncorrelated fit.
        verbose: Controls the amount of information written to the screen.

    Returns:
        The function returns the fitting parameters, the chi^2 and the p-value
        for every bootstrap sample.
    """
    try:
#        errfunc = lambda p, x, y, error: np.dot(error, (y-fitfunc(p,x)).T) + constrain(p)
        errfunc = lambda p, x, y, error: np.dot(error, (y-fitfunc(p,x)).T)
    
        # compute inverse, cholesky decomposed covariance matrix
        if not correlated:
            cov = np.diag(np.diagonal(np.cov(Y.T)))
        else:
            cov = np.cov(Y.T)
#        cov_inv = np.linalg.inv(cov)
#        print(cov_inv)
        cov = (np.linalg.cholesky(np.linalg.inv(cov))).T
        # degrees of freedom
        dof = float(Y.shape[1]-len(start_parm)) 
        # create results arrays
        res = np.zeros((Y.shape[0], len(start_parm)))
        res_cov = np.zeros((len(start_parm), len(start_parm)))
        chisquare = np.zeros(Y.shape[0])
        # The FIT to the boostrap samples
        for b in range(0, Y.shape[0]):
#        for b in range(1):
            p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm, args=(X, Y[b,:], cov), full_output=1, factor=0.1)
#            print b, mesg, ier, p
            chisquare[b] = float(sum(infodict['fvec']**2.))
            res[b] = np.array(p)
#            if b==0:
#               res_cov = cov1*chisquare[b]/dof
        # p-value calculated
        pvals_originfit = 1. - scipy.stats.chi2.cdf(chisquare[0], dof)
#        print(res.shape)
#        res_mean, res_std = af.calc_error(res)
#        print(res_mean)
#        print(res_std)
#        print(res[0])
#        print("chi2/dof = %.6f/%i = %.6f" % (chisquare[0], dof, chisquare[0]/dof))
#        print("pval = %.6f" % (pvals_originfit))
#        if (res[0,1] < 0.0) or (res[0,1] > 1.5):
#           return 0
#        else:
        return res, res_cov, chisquare[0], dof, pvals_originfit
    
    except:
        return 0

def fitting1(fitfunc, X1, X2, Y, start_parm, correlated=True, verbose=False):
    """A function that fits a correlation function.

    This function fits the given function fitfunc to the data given in X and Y.
    The function needs some start values, given in start_parm, and can use a
    correlated or an uncorrelated fit.

    Args:
        fitfunc: The function to fit to the data.
        X: The time slices.
        Y: The bootstrap samples of the data.
        start_parm: The starting parameters for the fit.
        correlated: Flag to use a correlated or uncorrelated fit.
        verbose: Controls the amount of information written to the screen.

    Returns:
        The function returns the fitting parameters, the chi^2 and the p-value
        for every bootstrap sample.
    """
    try:
        errfunc = lambda p, x1, x2, y, error: np.dot(error, (y-fitfunc(p,x1,x2)).T)
    
        # compute inverse, cholesky decomposed covariance matrix
        if not correlated:
            cov = np.diag(np.diagonal(np.cov(Y.T)))
        else:
            cov = np.cov(Y.T)
        cov = (np.linalg.cholesky(np.linalg.inv(cov))).T
        # degrees of freedom
        dof = float(Y.shape[1]-len(start_parm)) 
        # create results arrays
        res = np.zeros((Y.shape[0], len(start_parm)))
        res_cov = np.zeros((len(start_parm), len(start_parm)))
        chisquare = np.zeros(Y.shape[0])
        # The FIT to the boostrap samples
        for b in range(0, Y.shape[0]):
            p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm, args=(X1, X2, Y[b,:], cov), full_output=1, factor=0.1)
            chisquare[b] = float(sum(infodict['fvec']**2.))
            res[b] = np.array(p)
            if b==0:
               res_cov = cov1*chisquare[b]/dof
        # p-value calculated
        pvals_originfit = 1. - scipy.stats.chi2.cdf(chisquare[0], dof)
#        print(res.shape)
#        res_mean, res_std = af.calc_error(res)
#        print(res_mean)
#        print(res_std)
#        print(res[0])
#        print("chi2/dof = %.6f/%i = %.6f" % (chisquare[0], dof, chisquare[0]/dof))
#        print("pval = %.6f" % (pvals_originfit))
#        if (res[0,1] < 0.0) or (res[0,1] > 1.5):
#           return 0
#        else:
        return res, res_cov, chisquare[0], dof, pvals_originfit
    
    except:
        return 0

def weighted_quantile(data, weights, quantile):
    """Compute the weighted quantile, where a fixed percentage of the sum of
    all weights lie below.

    Args:
        data: A numpy-array of the data points the quantile is taken from.
        weights: A numpy-array containing the weights for each point in data. 
              Must be of same shape and have same order as data.
        quantile: The percentage of weights to be below the quantile. 
              0.5 is the weighted median
    """
    ind_sorted = np.argsort(data)
    sorted_data = data[ind_sorted]
    sorted_weights = weights[ind_sorted]
    # Compute the auxiliary arrays
    Sn = np.cumsum(sorted_weights)
    Pn = (Sn-0.5*sorted_weights)/np.sum(sorted_weights)
    # Get the value of the weighted median
    interpolated_quant = np.interp(quantile, Pn, sorted_data)

    return interpolated_quant


def set_fit_interval(nt, t_low, t_up, intervalsize):
    """Initialize intervals to fit in with borders given for every principal
    correlator

    Args: 
        data: The lattice results to fit to. Necessary to obtain the number of
              gevp-eigenvalues.
        lolist: List of lower interval borders for every gevp-eigenvalue.
        uplist: List of upper interval borders for every gevp-eigenvalue.
        intervallsize: Minimal number of points to be contained in the 
                interval

    Returns:
        fit_intervals: list of tuples (lo, up) for every gevp-eigenvalue.
    """
    fit_intervals = []
    if t_up > nt or t_low < 0:
       print("upper bound of fit range out of the time extent.")
       os.sys.exit(-1) 
    for lo in range(t_low, t_up):
        for up in range(t_low+1, t_up+1):
            if (up - lo + 2) > intervalsize:
                fit_intervals.append((lo, up))

    return fit_intervals

def genfit(data, tmin, tmax, min_t_interval, fitfunc, start_params, out_originfit, out_bootstrapfit, out_finalresults, out_all, verbose=False):
    """Fit and plot the correlation function.
    
    Args:
        data: The correlation functions.
        fit_intervalls: List of intervalls for the fit for the different
              correlation functions.
        fitfunc: The function to fit to the data.
        start_params: The starting parameters for the fit function.
        verbose: Amount of information printed to screen.

    Returns:
        res: Result of the fit to each bootstrap sample.
        chi2: Chi^2 for every fit
        pval: p-value for every fit.
    """
    # init variables
    nboot = data.shape[0]
    Lt = data.shape[1]
    T2 = Lt/2
    npar = len(start_params)
    fit_intervals = np.asarray(set_fit_interval(T2, tmin, tmax, min_t_interval))
    ninter = len(fit_intervals)
#    print(fit_intervals)
    print(ninter)

    # set fit data
    tlist = np.linspace(0., float(T2), float(T2), endpoint=False)

    # initialize array for every principal correlator
    res = []  #ninter, nboot, npar
    res_cov = [] #ninter, npar*npar
    chi2 = []  #nboot
    pval = []  #nboot
    dof = []
    failed_inter = []
    weighted_res = np.zeros((nboot, npar))
    stat_error = np.zeros((npar))
    sys_error = np.zeros((npar, 2))
    for _i in range(ninter):
            print(_i)
            lo, up = fit_intervals[_i]
            if verbose:
                print("Interval [%d, %d]" % (lo, up))
                print("fitting correlation function")
                print(tlist[lo:up+1])
            if fitting(fitfunc, tlist[lo:up+1], data[:,lo:up+1], start_params, verbose=False) == 0:
               print("fit range [%i, %i] failed. continue..." % (lo, up))
               failed_inter.append(_i)
            else:
               _res, _res_cov, _chi2, _dof, _pval = fitting(fitfunc, tlist[lo:up+1], data[:,lo:up+1], start_params, verbose=False)
               res.append(_res)
               res_cov.append(_res_cov)
               chi2.append(_chi2)
               pval.append(_pval)
               dof.append(_dof)

    fit_intervals = np.delete(fit_intervals, failed_inter, 0)
    ninter = len(fit_intervals)
    print("n intervals after removing failed intervals:")
    print(ninter)
    res = np.array(res)
    res_cov = np.array(res_cov)
    chi2 = np.array(chi2)
    pval = np.array(pval)

    res_error = np.zeros((ninter, npar))
    for _i in range(ninter):
       for _j in range(npar):
          res_error[_i, _j] = np.std(res[_i, :, _j])
    
    weight = np.zeros((ninter, npar))
    for _i in range(ninter):
        for _j in range(npar):
           weight[_i, _j] = ((1.-2.*np.fabs(pval[_i] - 0.5))*np.amin(res_error[:,_j])/res_error[_i,_j])**2

    for _i in range(nboot):
      for _j in range(npar):
#	weighted_res[_i, _j] = weighted_quantile(res[:,_i, _j], weight[:,_j], 0.5) 
	weighted_res[_i, _j] = weighted_quantile(res[:,_i, _j], weight[:,1], 0.5) 

    for _j in range(npar):
	stat_error[_j] = np.std(weighted_res[:, _j])
	sys_error[_j, 0] = weighted_res[0, _j] - weighted_quantile(res[:, 0, _j], weight[:, 1], 0.16)
	sys_error[_j, 1] = weighted_quantile(res[:, 0, _j], weight[:, 1], 0.84) - weighted_res[0, _j]
#	sys_error[_j, 0] = weighted_res[0, _j] - weighted_quantile(res[:, 0, _j], weight[:, _j], 0.16)
#	sys_error[_j, 1] = weighted_quantile(res[:, 0, _j], weight[:, _j], 0.84) - weighted_res[0, _j]

#    out_origin_fit = np.column_stack((fit_intervals, res[:,0,:], res_cov, chi2, pval, weight))

#    np.savetxt(out_originfit, out_origin_fit, fmt=["%i", "%i", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f"])
    np.savetxt(out_bootstrapfit, weighted_res, fmt=["%.16f", "%.16f"] )

# output of original fit for plots
    np.savez(out_originfit, fitranges=fit_intervals, dof=dof, par=res[:,0,:], cov=res_cov, chi2=chi2, pvals=pval, weight=weight) 

    out_final = np.zeros((npar, 5))
    for i in range(npar):
       out_final[i] = np.array([weighted_res[0, i], np.mean(weighted_res[:,i]), stat_error[i], sys_error[i, 0], sys_error[i, 1]])
    np.savetxt(out_finalresults, out_final, fmt=["%.16f", "%.16f", "%.16f", "%.16f", "%.16f"])

    for _i in range(npar):
      np.savez("%s_par%i.npz" % (out_all, _i), weight=weight[:,_i], par=res[:,:,_i])

