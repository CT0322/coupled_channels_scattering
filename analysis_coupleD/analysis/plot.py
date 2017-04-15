
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import analyze_fcts as af

def plot_data_with_fit(X, Y, dY, fitfunc, args, plotrange, T, fitranges, pdfplot, pval, chi2, dof, weight, logscale=False, xlim=None, ylim=None, addpars=False):
    """A function that plots data and the fit to the data.

    The plot is saved to pdfplot. It is assumed that pdfplot is a pdf backend to
    matplotlib so that multiple plots can be saved to the object.

    Args:
        X: The data for the x axis.
        Y: The data for the y axis.
        dY: The error on the y axis data.
        fitfunc: The function to fit to the data.
        args: The parameters of the fit function from the fit.
        plotrange: A list with two entries, the lower and upper range of the
                   plot.
        pdfplot: A PdfPages object in which to save the plot.
        logscale: Make the y-scale a logscale.
        xlim, ylim: limits for the x and y axis, respectively
        setLimits: Set limits to the y range of the plot.
        fitrange: A list with two entries, bounds of the fitted function.
        addpars: if there are additional parameters for the fitfunction 
                 contained in args, set to true
        pval: write the p-value in the plot if given

    Returns:
        Nothing.
    """
    # check boundaries for the plot
    if isinstance(plotrange, (np.ndarray, list, tuple)):
        plotrange = np.asarray(plotrange).flatten()
        if plotrange.size < 2:
            raise IndexError("plotrange has not enough indices")
        else:
            l = int(plotrange[0])
            u = int(plotrange[1])
    else:
        l = 0
        u = T/2 + 1
    ninter = fitranges.shape[0]
    p1 = plt.errorbar(X[l:u+1], Y[l:u+1], dY[l:u+1], fmt='x' + 'b')
    for _n in range(ninter):   
    # plotting the fit function, check for seperate range
        lfunc = int(fitranges[_n,0])
        ufunc = int(fitranges[_n,1])
        x1 = np.arange(l, u, 0.01)
        lfunci = np.where(x1==lfunc)[0]
        ufunci = np.where(x1==ufunc)[0]
        ui = x1.shape[0]
        y1 = []
        if addpars:
            for i in x1:
            # the star in front of the args is needed
                y1.append(fitfunc(args[_n, 0],i,*args[_n, 1:]))
        else:    
            for i in x1:
                y1.append(fitfunc(args[_n],i))
        y1 = np.asarray(y1)

        p2, = plt.plot(x1[lfunci:ufunci], y1[lfunci:ufunci], "r")
        if l<lfunc:
            p3, = plt.plot(x1[0:lfunci], y1[0:lfunci], "g")
        if ufunc<u:
            p4, = plt.plot(x1[ufunci:ui], y1[ufunci:ui], "g")
        # adjusting the plot style
        plt.title("$\chi^{2}$=%.2f, p=%.6f, w=%.6f, m=%.6f" % (chi2[_n]/dof[_n], pval[_n], weight[_n, 1], args[_n,1]))
        plt.xlabel('t')
        plt.ylabel('C(t)')
        if logscale:
            plt.yscale("log")
        # set the axis ranges
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        # save pdf
        pdfplot.savefig()
        p2.remove()
        if l<lfunc:
            p3.remove()
        if ufunc<u:
            p4.remove()
    plt.clf()
    pdfplot.close()

def plot_data_with_fit_mean(X, Y, dY, fitfunc, fitrange, args, masserror, plotrange, T, pdfplot, logscale=False, xlim=None, ylim=None, addpars=False):
    """A function that plots data and the fit to the data.

    The plot is saved to pdfplot. It is assumed that pdfplot is a pdf backend to
    matplotlib so that multiple plots can be saved to the object.

    Args:
        X: The data for the x axis.
        Y: The data for the y axis.
        dY: The error on the y axis data.
        fitfunc: The function to fit to the data.
        args: The parameters of the fit function from the fit.
        plotrange: A list with two entries, the lower and upper range of the
                   plot.
        pdfplot: A PdfPages object in which to save the plot.
        logscale: Make the y-scale a logscale.
        xlim, ylim: limits for the x and y axis, respectively
        setLimits: Set limits to the y range of the plot.
        fitrange: A list with two entries, bounds of the fitted function.
        addpars: if there are additional parameters for the fitfunction 
                 contained in args, set to true
        pval: write the p-value in the plot if given

    Returns:
        Nothing.
    """
    # check boundaries for the plot
    if isinstance(plotrange, (np.ndarray, list, tuple)):
        plotrange = np.asarray(plotrange).flatten()
        if plotrange.size < 2:
            raise IndexError("plotrange has not enough indices")
        else:
            l = int(plotrange[0])
            u = int(plotrange[1])
    else:
        l = 0
        u = T/2 + 1
    fitl = fitrange[0]
    fitu = fitrange[1]
    p1 = plt.errorbar(X[l:u+1], Y[l:u+1], dY[l:u+1], fmt='x' + 'b')
    x1 = np.arange(l, u+0.01, 0.01)
    lfunci = np.where(x1==l)[0]
    ufunci = np.where(x1==u)[0]
    fitlfunci = np.where(x1==fitl)[0]
    fitufunci = np.where(x1==fitu)[0]
#    print(lfunci)
#    print(ufunci)
    y1 = []
    for i in x1:
        y1.append(fitfunc(args,i))
    y1 = np.asarray(y1)

    p2 = plt.plot(x1[lfunci:fitlfunci], y1[lfunci:fitlfunci], "k")
    p3 = plt.plot(x1[fitlfunci:fitufunci], y1[fitlfunci:fitufunci], "r")
    p4 = plt.plot(x1[fitufunci:ufunci], y1[fitufunci:ufunci], "k")
        # adjusting the plot style
    plt.xlabel('t')
    plt.ylabel('C(t)')
    plt.title("m = %.6f ($\pm$ %.6f) (- %.6f / + %.6f)" % (args[1], masserror[0], masserror[1], masserror[2]))
    if logscale:
        plt.yscale("log")
        # set the axis ranges
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
        # save pdf
    pdfplot.savefig()
    plt.clf()
    pdfplot.close()



def plot_data_with_fit_1(X, Y, dY, fitfunc, fitrange, args, masserror, chi2, dof, pval, plotrange, T, pdfplot, logscale=False, xlim=None, ylim=None, addpars=False):
    """A function that plots data and the fit to the data.

    The plot is saved to pdfplot. It is assumed that pdfplot is a pdf backend to
    matplotlib so that multiple plots can be saved to the object.

    Args:
        X: The data for the x axis.
        Y: The data for the y axis.
        dY: The error on the y axis data.
        fitfunc: The function to fit to the data.
        args: The parameters of the fit function from the fit.
        plotrange: A list with two entries, the lower and upper range of the
                   plot.
        pdfplot: A PdfPages object in which to save the plot.
        logscale: Make the y-scale a logscale.
        xlim, ylim: limits for the x and y axis, respectively
        setLimits: Set limits to the y range of the plot.
        fitrange: A list with two entries, bounds of the fitted function.
        addpars: if there are additional parameters for the fitfunction 
                 contained in args, set to true
        pval: write the p-value in the plot if given

    Returns:
        Nothing.
    """
    # check boundaries for the plot
    if isinstance(plotrange, (np.ndarray, list, tuple)):
        plotrange = np.asarray(plotrange).flatten()
        if plotrange.size < 2:
            raise IndexError("plotrange has not enough indices")
        else:
            l = int(plotrange[0])
            u = int(plotrange[1])
    else:
        l = 0
        u = T/2 + 1
    fitl = fitrange[0]
    fitu = fitrange[1]
    p1 = plt.errorbar(X[l:u+1], Y[l:u+1], dY[l:u+1], fmt='x' + 'b')
    x1 = np.arange(l, u+0.01, 0.01)
    lfunci = np.where(x1==l)[0]
    ufunci = np.where(x1==u)[0]
    fitlfunci = np.where(x1==fitl)[0]
    fitufunci = np.where(x1==fitu)[0]
#    print(lfunci)
#    print(ufunci)
    y1 = []
    for i in x1:
        y1.append(fitfunc(args,i))
    y1 = np.asarray(y1)

    p2 = plt.plot(x1[lfunci:fitlfunci], y1[lfunci:fitlfunci], "k")
    p3 = plt.plot(x1[fitlfunci:fitufunci], y1[fitlfunci:fitufunci], "r")
    p4 = plt.plot(x1[fitufunci:ufunci], y1[fitufunci:ufunci], "k")
        # adjusting the plot style
    plt.xlabel('t')
    plt.ylabel('C(t)')
#    plt.title("m = %.6f (%.6f), $\chi^{2}$/dof=%.2f/%i, p=%.6f " % (args[0], masserror, chi2, dof, pval))
    plt.title("m = %.6f (%.6f), $\chi^{2}$/dof=%.2f/%i, p=%.6f " % (args[1], masserror, chi2, dof, pval))
    if logscale:
        plt.yscale("log")
        # set the axis ranges
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
        # save pdf
    pdfplot.savefig()
    plt.clf()
    pdfplot.close()

# this can be used to plot the chisquare distribution of the fits
#  x = np.linspace(scipy.stats.chi2.ppf(1e-6, dof), scipy.stats.chi2.ppf(1.-1e-6, dof), 1000)
#  hist, bins = np.histogram(chisquare, 50, density=True)
#  width = 0.7 * (bins[1] - bins[0])
#  center = (bins[:-1] + bins[1:]) / 2
#  plt.xlabel('x')
#  plt.ylabel('chi^2(x)')
#  plt.grid(True)
#  plt.plot(x, scipy.stats.chi2.pdf(x, dof), 'r-', lw=2, alpha=1, label='chi2 pdf')
#  plt.bar(center, hist, align='center', width=width)
#  plt.show()

def plot_data(X, Y, dY, pdfplot, plotrange=None, logscale=True, xlim=None, ylim=None):
    """A function that plots a correlation function.

    This function plots the given data points and the fit to the data. The plot
    is saved to pdfplot. It is assumed that pdfplot is a pdf backend to
    matplotlib so that multiple plots can be saved to the object.

    Args:
        X: The data for the x axis.
        Y: The data for the y axis.
        dY: The error on the y axis data.
        pdfplot: A PdfPages object in which to save the plot.
        plotrange: A list with two entries, the lower and upper range of the
                   plot.
        logscale: Make the y-scale a logscale.
        xlim: tuple of the limits on the x axis
        ylim: tuple of the limits on the y axis

    Returns:
        Nothing.
    """
    # check boundaries for the plot
    if isinstance(plotrange, (np.ndarray, list, tuple)):
        plotrange = np.asarray(plotrange).flatten()
        if plotrange.size < 2:
            raise indexerror("plotrange is too small")
        else:
            l = int(plotrange[0])
            u = int(plotrange[1])
    else:
        l = 0
        u = x.shape[0]

    # plot the data
    p1 = plt.errorbar(X[l:u], Y[l:u], dY[l:u], fmt='x' + 'b', label="data")

    # adjusting the plot style
    plt.grid(True)
    #plt.xlabel(label[1])
    #plt.ylabel(label[2])
    #plt.title(label[0])
    plt.legend()
    if logscale:
        plt.yscale('log')
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)

    # save pdf and clear plot
    pdfplot.savefig()
    plt.clf()

    return

#def plot_histogram(data, data_weight, lattice, d, label, path="./plots/", 
#                   plotlabel="hist", verbose=True):
def plot_histogram(data, data_weight, histplot):
    """Plots histograms for the given data set.

    The function plots the weighted distribution of the data, the unweighted
    distribution and a plot containing both the weighted and the unweighted
    distribution.

    Args:
        data: Numpy-array of fit values for mulitple fit intervalls. Will be 
              depicted on x-axis.
        data_weight: The weights corresponding to data. Must have same shape
              and order as data. Their sum per bin is the bin height.
        lattice: The name of the lattice, used for the output file.
        d:    The total momentum of the reaction.
        label: Labels for the title and the axis.
        path: Path to the saving place of the plot.
        plotlabel: Label for the plot file.
        verbose: Amount of information printed to screen.

    Returns:
    """
    ninter = data.shape[0]


    # The histogram
    # generate weighted histogram
    hist, bins = np.histogram(data, 50, weights=data_weight, density=True)
    # generate the unweighted histogram
    uhist, ubins = np.histogram(data, 50, weights=np.ones_like(data_weight),
                                density=True)

    # prepare the plot for the weighted histogram
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    # set labels for axis
#    plt.ylabel('weighted distribution of ' + label[2])
#    plt.title('fit methods individually')
#    plt.grid(True)
    # plot
    plt.bar(center, hist, align='center', width=width, color='r', alpha=0.5,
            label='weighted data')
    # save and clear
    histplot.savefig()
    plt.clf()

    # prepare plot for unweighted histogram
    # the center and width stays the same for comparison
#    plt.ylabel('unweighted distribution of ' + label[2])
#    plt.title('fit methods individually')
#    plt.grid(True)
    # plot
    plt.bar(center, uhist, align='center', width=width, color='b', alpha=0.5,
            label='unweighted data')

    # save and clear
    histplot.savefig()
    plt.clf()

    # plot both histograms in same plot
#    plt.ylabel('distribution of ' + label[2])
#    plt.title('fit methods individually')
#    plt.grid(True)
    # plot
    plt.bar(center, hist, align='center', width=width, color='r', alpha=0.5,
            label='weighted data')
    plt.bar(center, uhist, align='center', width=width, color='b', alpha=0.5,
            label='unweighted data')
    plt.legend()

    # save and clear
    histplot.savefig()
    plt.clf()

    # close plotfile
    histplot.close()
    return

