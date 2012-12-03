#!/usr/bin/env python
from __future__ import division
import itertools
import pylab
from numpy import *
import lomb
import math
from operator import itemgetter, mul
#from bintrees import FastBinaryTree as Tree
from six.moves import filter, zip

DATA_FILES = {
    "IGR_J18027-2016_17-30": "log_IGRJ18027-2016_17-30.dat",
    "GX13+1_30-60": "log_GX13+1_30-60.dat"
}


def quadrature(dx):
    sqerrs = map(lambda dx: dx ** 2, dx)
    return math.sqrt(math.fsum(sqerrs))


def fractional(x, dx):
    sqfracerrs = map(lambda x, dx: (dx / x) ** 2, x, dx)
    mult = reduce(mul, x)
    return mult * math.sqrt(math.fsum(sqfracerrs))


def weighted_av(x, dx):
    weights = map(lambda x: 1 / (x ** 2), dx)
    xw = map(lambda x, w: x * w, x, weights)
    return math.fsum(xw) / math.fsum(weights)


def weighted_av_error(x, dx):
    def fractxy(x, dx, y):
        return y * (math.sqrt((dx / x) ** 2))
    w = map(lambda x: 1 / (x ** 2), dx)
    q = map(lambda x, w: x * w, x, w)
    dq = map(fractxy, x, dx, q)
    return fractxy(math.fsum(q), quadrature(dq), weighted_av(x, dx))


def rebin_error(rows, keys=("flux", "flux_error")):
    flux = []
    error = []
    for row in rows:
        flux.append(row[keys[0]])
        error.append(row[keys[1]])
    return weighted_av(flux, error), weighted_av_error(flux, error)


def import_data(file_name):
    """
    Imports text file of format specified by T.Bird, outputs into arrays:
    MJD, flux, flux_error, OOA, exposure.
    example MJD, flux, flux_error, OOA, exposure = import_data(file_name)
    """
    return loadtxt(file_name, usecols=(1, 2, 3, 5, 6), dtype={
        "names": ("MJD", "flux", "flux_error", "OAA", "exposure"),
        "formats": ("f", "f", "f", "f", "f")
    })


def rebin(rows, bin_size=30):
    for bin_key, bin in itertools.groupby(rows, lambda row: (row["MJD"] // bin_size) * bin_size):
        rows = tuple(bin)
        flux, error = rebin_error(rows)
        yield (bin_key, flux, error)


def plot_rebin(data):
    x, y, errors = zip(*rebin(data))
    x_shift = zeros(len(x))
    for i in range(len(x)):
        x_shift[i] = x[i] - x[0]
    pylab.errorbar(x_shift, y, yerr=errors, fmt=".")
    pylab.ylim([-20,20])
    pylab.xlabel("time")
    pylab.ylabel("counts")
    pylab.title('IGRJ18027 2016 17-30 KeV 1 Month Lightcurve')
    pylab.show()


def lombscargle(rows, oversample=6., nyquist=10, keys=("MJD", "flux")):
    x = []
    y = []
    for row in rows:
        x.append(row[keys[0]])
        y.append(row[keys[1]])

    fx, fy, nout, jmax, prob = lomb.fasper(array(x), array(y), oversample, nyquist)

    return (zip(fx, fy), jmax)


def periodogram(data, plot=False):
    data, jmax = lombscargle(data)

    # calculate the reciprocal of the x dimension
    data = list(itertools.starmap(lambda x, y: (1 / x, y), data))
    period, power = zip(*data)
    period_max = data[jmax][0]

    print period_max

    if plot == True:
        pylab.vlines(period, 0, array(power), color='k', linestyles='solid')
        pylab.xlabel("period, days")
        pylab.ylabel("power")
        pylab.xlim([0, 40])
        pylab.title('IGRJ18027 2016_17-30 KeV Periodogram ')
        pylab.show()

    return period_max


def fold(rows, period, number_bins=10):
    for row in rows:
        row["MJD"] %= period
    rows.sort(key=itemgetter("MJD"))

    MJD_bins = {k: list() for k in linspace(0, period_max, number_bins)[:-1]}  # dropping the last element

    for row in rows:
        # discover the biggest key that is smaller than row["MJD"]
        biggest_key = -1
        for key in MJD_bins.keys():
            if (key > biggest_key) and (key <= row["MJD"]):
                biggest_key = key

        MJD_bins[biggest_key].append(row)

    # weighted average of fluxes in bins
    for key, value in MJD_bins.items():
        if value:
            flux, error = rebin_error(value)
            yield (key, flux, error)

def plot_fold(data, period_max):
    #crap = list(fold(data, period_max))
    x, y, errors = zip(*fold(data, period_max, 30))
    x = itertools.chain.from_iterable([x, map(lambda x: x + period_max, x)])
    y = itertools.chain.from_iterable([y, y])
    errors = itertools.chain.from_iterable([errors, errors])
    #pylab.plot(list(x), list(y), '.')
    pylab.errorbar(list(x), list(y), list(errors), fmt='.')
    pylab.xlabel("time")
    pylab.ylabel("counts")
    pylab.title('IGRJ18027 2016 17-30 KeV Folded Lightcurve')
    pylab.show()


if __name__ == "__main__":
    data = list(filter(
        lambda row: (row["OAA"] < 12) and (row["exposure"] > 500),
        import_data(DATA_FILES["IGR_J18027-2016_17-30"])
    ))

    plot_rebin(data)
    period_max = periodogram(data)
    plot_fold(data, period_max)
