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

def fitting(fitfunc, X, Y, start_parm, correlated=True, verbose=True):
    """A function that fits a correlation function.

    This function fits the given function fitfunc to the data given in X and Y.
    The function needs some start values, given in start_parm, and can use a
    correlated or an uncorrelated fit.

    Args:
        fitfunc: The function to fit to the data.
        X: The time slices.
        Y: The bootstrap samples of the data.
        start_parm: The starting parameters for the fit.
        E_single: single particle energies entering the ratio R
        correlated: Flag to use a correlated or uncorrelated fit.
        verbose: Controls the amount of information written to the screen.

    Returns:
        The function returns the fitting parameters, the chi^2 and the p-value
        for every bootstrap sample.
    """
    errfunc = lambda p, x, y, error: np.dot(error, (y-fitfunc(p,x)).T)
    
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
        p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm, 
            args=(X, Y[b,:], cov), full_output=1, factor=0.1)
        chisquare[b] = float(sum(infodict['fvec']**2.))
        res[b] = np.array(p)
        if b==0:
    #       print(cov1)
           res_cov = cov1*chisquare[b]/dof
    #       print(res_cov)	
    # calculate mean and standard deviation
    res_mean, res_std = af.calc_error(res)
    # chi2 = np.median(chisquare)
    # p-value calculated
    pvals_originfit = 1. - scipy.stats.chi2.cdf(chisquare[0], dof)
    
    # The fit to the mean value
    y = np.mean(Y, axis=0)
    p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm, \
                                   args=(X, y, cov), full_output=1)
    chisquare_meanfit = float(sum(infodict['fvec']**2.))
    pvals_meanfit = 1. - scipy.stats.chi2.cdf(chisquare_meanfit, dof)
    # writing results to screen
    if verbose:
        if correlated:
            print("fit results for a correlated fit:")
        else:
            print("fit results for an uncorrelated fit:")
        print("degrees of freedom: %f\n" % dof)
        
        print("bootstrap fit:")
        for rm, rs in zip(res_mean, res_std):
            print("  %.6e +/- %.6e" % (rm, rs))
        #print("Chi^2/dof: %.6e +/- %.6e\n" % (chi2/dof, np.std(chisquare)/dof))

        print("mean value fit:")
        for rm, rs in zip(p, res_std):
            print("  %.6e +/- %.6e" % (rm, rs))
        print("  Chi^2/dof: %.6e " % (chisquare_meanfit / dof))
        print("  p-value: %lf" % pvals_meanfit) 

        print("original data fit:")
        for rm, rs in zip(res[0], res_std):
            print("  %.6e +/- %.6e" % (rm, rs))
        print("  Chi^2/dof: %.6e " % (chisquare[0]/dof))
        print("  p-value: %lf" % pvals_originfit) 
    return res, res_cov.flatten(), chisquare[0]/dof, pvals_originfit
#    return res, chisquare[0]/dof, pvals_originfit

#def quantile_1D(data, weights, quantile):
#    ind_sort = np.argsort(data)
#    sort_data = data[ind_sort]
#    sort_weig = wheights[ind_sort]
#    Sn = np.cumsum(sort_weig)
#    Pn = (Sn - 0.5*sort_weig) / np.sum(sort_weig)
#    return np.interp(quantile, Pn, sort_data)


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
    if t_up > nt - 1:
       print("upper bound for fit greater than time extent of data")
       print("using data time extend!")
       t_up = nt - 1
    for lo in range(t_low, t_up):
        for up in range(t_low+1, t_up+1):
            if (up - lo + 2) > intervalsize:
                fit_intervals.append((lo, up))

    return fit_intervals

def genfit(data, tmin, tmax, min_t_interval, fitfunc, start_params, out_originfit, out_bootstrapfit, out_finalresults, verbose=True):
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
    print(fit_intervals)
    print(ninter)

    # set fit data
    tlist = np.linspace(0., float(T2), float(T2), endpoint=False)

    # initialize array for every principal correlator
    res = np.zeros((nboot, npar, ninter))
    res_error = np.zeros((npar, ninter))
    res_cov = np.zeros((ninter, npar*npar))
    chi2 = np.zeros((ninter))
    pval = np.zeros((ninter))
    weight = np.zeros((npar, ninter))
    weighted_res = np.zeros((nboot, npar))
    stat_error = np.zeros((npar))
    sys_error = np.zeros((npar, 2))

    for _i in range(ninter):
            lo, up = fit_intervals[_i]
            if verbose:
                print("Interval [%d, %d]" % (lo, up))
                print("fitting correlation function")
                print(tlist[lo:up+1])
            res[:,:,_i], res_cov[_i,:], chi2[_i], pval[_i] = fitting(fitfunc, tlist[lo:up+1], data[:,lo:up+1], start_params, verbose=True)
            if verbose:
                print("p-value: %.7lf\nChi^2/dof: %.7lf\nresults:"
                      % (pval[_i], chi2[_i]/( (up - lo + 1) - len(start_params))))
                for p in enumerate(res[0,:,_i]):
                    print("\tpar %d = %lf" % p)
                print(" ")

    for _i in range(ninter):
       for _j in range(npar):
          res_error[_j, _i] = np.std(res[:,_j,_i])
    
    for _i in range(ninter):
        for _j in range(npar):
           weight[_j, _i] = (1.-2.*np.fabs(pval[_i] - 0.5))*np.amax(res_error[_j,:])/res_error[_j,_i]

    for _i in range(nboot):
      for _j in range(npar):
	weighted_res[_i, _j] = weighted_quantile(res[_i, _j, :], weight[_j, :], 0.5) 

    for _j in range(npar):
	stat_error[_j] = np.std(weighted_res[:, _j])
	sys_error[_j, 0] = weighted_res[0, _j] - weighted_quantile(res[0, _j, :], weight[_j, :], 0.16)
	sys_error[_j, 1] = weighted_quantile(res[0, _j, :], weight[_j, :], 0.84) - weighted_res[0, _j]

    out_origin_fit = np.column_stack((fit_intervals, res[0,:,:].T, res_cov, chi2, pval, weight.T))
    np.savetxt(out_originfit, out_origin_fit, fmt=["%i", "%i", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f", "%.10f"])
    np.savetxt(out_bootstrapfit, weighted_res, fmt=["%.10f", "%.10f"] )

    out_final = np.zeros((npar, 4))
    for i in range(npar):
       out_final[i] = np.array([weighted_res[0, i], stat_error[i], sys_error[i, 0], sys_error[i, 1]])
    np.savetxt(out_finalresults, out_final, fmt=["%.10f", "%.10f", "%.10f", "%.10f"])
#    print("********************** fitted parameters: ****************************")
#    print(res[0])
#    print("********************** weight: ****************************")
#    print(weight)
#    print("********************** weighted parameters: ****************************")
#    print(weighted_res)
#    print("********************** stat and sys error: ****************************")
#    print(stat_error)
#    print(sys_error)
#    return res, chi2, pval, weight, weighted_res
#    return res[0], chi2, pval, weight, weighted_res, stat_error, sys_error  


