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
from .input_output import read_fitresults, write_fitresults

def fitting(fitfunc, X, Y, start_parm, E_single=None, correlated=True, verbose=True):
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
    if E_single is None:
        errfunc = lambda p, x, y, error: np.dot(error, (y-fitfunc(p,x)).T)
    else:
        errfunc = lambda p, x, y, e, error: np.dot(error, (y-fitfunc(p,x,e)).T)
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
    chisquare = np.zeros(Y.shape[0])
    # The FIT to the boostrap samples
    if E_single is None:
        for b in range(0, Y.shape[0]):
            p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm, 
                args=(X, Y[b,:], cov), full_output=1, factor=0.1)
            chisquare[b] = float(sum(infodict['fvec']**2.))
            res[b] = np.array(p)
    else:
        for b in range(0, Y.shape[0]):
            p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm,
                                       args=(X, Y[b,:],E_single[b], cov ), 
                                       full_output=1, factor=0.1)
            chisquare[b] = float(sum(infodict['fvec']**2.))
            res[b] = np.array(p)
    # calculate mean and standard deviation
    res_mean, res_std = af.calc_error(res)
    chi2 = np.median(chisquare)
    # p-value calculated
    pvals = 1. - scipy.stats.chi2.cdf(chisquare, dof)

    # The fit to the mean value
    y = np.mean(Y, axis=0)
    if E_single is None:
        p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm, \
                                   args=(X, y, cov), full_output=1)
    else:
        e_single = np.mean(E_single)
        p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm, \
                                   args=(X, y, e_single, cov), full_output=1)
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
        print("Chi^2/dof: %.6e +/- %.6e\n" % (chi2/dof,
              np.std(chisquare)/dof))

        print("mean value fit:")
        for rm, rs in zip(p, res_std):
            print("  %.6e +/- %.6e" % (rm, rs))
        print("Chi^2/dof: %.6e +/- %.6e\n" % (float(sum(infodict['fvec']**2.) /
              dof), np.std(chisquare)/dof))

        print("original data fit:")
        for rm, rs in zip(res[0], res_std):
            print("  %.6e +/- %.6e" % (rm, rs))
        print("Chi^2/dof: %.6e +/- %.6e" % (chisquare[0]/dof, np.std(chisquare)
              /dof))
        print("p-value: %lf" % pvals[0]) 

#    print("res:", res[0:10])
#    print("chisquare:", chisquare[0:10])
#    print("pvals:", pvals[0:10])
    print("test print pvals:", pvals)
    return res, chisquare, pvals

def quantile_1D(data, weights, quantile):
    ind_sort = np.argsort(data)
    sort_data = data[ind_sort]
    sort_weig = wheights[ind_sort]
    Sn = np.cumsum(sort_weig)
    Pn = (Sn - 0.5*sort_weig) / np.sum(sort_weig)
    return np.interp(quantile, Pn, sort_data)


def set_fit_interval(_data, lolist, uplist, intervalsize):
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
    data = np.atleast_3d(_data)
    ncorr = data.shape[2]
    fit_intervals = []
    for _l in range(ncorr):
        fit_intervals.append([])
        if uplist[_l] > data.shape[1] - 1:
            print("upper bound for fit greater than time extent of data")
            print("using data time extend!")
            uplist[_l] = data.shape[1] - 1
        for lo in range(lolist[_l], uplist[_l] + 1):
            for up in range(lolist[_l], uplist[_l] + 1):
                # the +2 comes from the fact that the interval contains one
                # more number than the difference between the boundaries and
                # because python would exclude the upper boundary but we
                # include it explicitly
                if (up - lo + 2) > intervalsize:
                    fit_intervals[_l].append((lo, up))

    return fit_intervals


def genfit(_data, fit_intervals, fitfunc, start_params, verbose=True):
    """Fit and plot the correlation function.
    
    Args:
        _data: The correlation functions.
        fit_intervalls: List of intervalls for the fit for the different
              correlation functions.
        fitfunc: The function to fit to the data.
        start_params: The starting parameters for the fit function.
        olddata: if not None, reuses old data at the location specified, if possible
        verbose: Amount of information printed to screen.

    Returns:
        res: Result of the fit to each bootstrap sample.
        chi2: Chi^2 for every fit
        pval: p-value for every fit.
    """
    data = np.atleast_3d(_data)
    # init variables
    nboot = data.shape[0]
    T2 = data.shape[1]
    ncorr = data.shape[2]
    npar = len(start_params)
    ninter = [len(fitint) for fitint in fit_intervals]
    # set fit data
    tlist = np.linspace(0., float(T2), float(T2), endpoint=False)
    # initialize empty arrays
    res = []
    chi2 = []
    pval = []
    # initialize array for every principal correlator
    for _l in range(ncorr):
        res.append(np.zeros((nboot, npar, ninter[_l])))
        chi2.append(np.zeros((nboot, ninter[_l])))
        pval.append(np.zeros((nboot, ninter[_l])))
    for _l in range(ncorr):
        # setup
        mdata, ddata = af.calc_error(data[:,:,_l])
        for _i in range(ninter[_l]):
            lo, up = fit_intervals[_l][_i]
            if verbose:
                print("Interval [%d, %d]" % (lo, up))
                print("correlator %d" % _l)
                print("fitting correlation function")
                print(tlist[lo:up+1])
            res[_l][:,:,_i], chi2[_l][:,_i], pval[_l][:,_i] = fitting(fitfunc, tlist[lo:up+1], data[:,lo:up+1,_l], start_params, verbose=False)
            if verbose:
                print("p-value %.7lf\nChi^2/dof %.7lf\nresults:"
                      % (pval[_l][ 0, _i], chi2[_l][0,_i]/( (up - lo + 1) - len(start_params))))
                for p in enumerate(res[_l][0,:,_i]):
                    print("\tpar %d = %lf" % p)
                print(" ")
    return res, chi2, pval


def compute_weight(corr, params):
    """compute the weights for the histogram

    Args:
        corr: the correlation function
        params: the p-values for each correlator

    Returns:
        list of weigths if len(params)!=0, else empty list
    """
    errors = np.std(corr, axis=1)
    max_err = np.amax(errors)
    weights = []
    if len(params) != 0:
        for i in range(0, params.shape[0]):
            w = (1.-2*abs(params[i,1]-0.5))*max_err/errors[i]
            weights.append(w**2)
    return weigths

