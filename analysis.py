#!/usr/bin/env python
from __future__ import division
import itertools
import pylab
from numpy import *
import lomb
import math
from operator import itemgetter, mul
from bintrees import FastBinaryTree as Tree
from six.moves import filter, zip
import random
from multiprocessing import Pool
from collections import defaultdict

DATA_FILES = {
    "IGR J18027-2016 17-30": "log_IGRJ18027-2016_17-30.dat",
    "IGR J18027-2016 30-60": "log_IGRJ18027-2016_30-60.dat",
    "GX13+1_30-60": "log_GX13+1_30-60.dat",
    "1A0535+262 18-60": "log_1A0535+262_18-60.dat",
    "AXJ1910.7+0917 18-60": "log_AXJ1910.7+0917_18-60.dat",
    "AXJ1910.7+0917 100-300": "log_AXJ1910.7+0917_100-300.dat",
    "CenX-3 18-60": "log_CenX-3_18-60.dat",
    "EXO2030+375 18-60": "log_EXO2030+375_18-60.dat",
    "IGRJ08408-4503 18-60": "log_IGRJ08408-4503_18-60.dat",
    "SMCX-1 18-60": "log_SMCX-1_18-60.dat",
    "OAO1657-415 18-60": "log_OAO1657-415_18-60.dat",
    "IGR J11215-5952 18-60": "log_IGRJ11215-5952_18-60.dat",
    "IGR J16418-4532 18-60": "log_IGRJ16418-4532_18-60.dat",
    "IGR J16465-4507 18-60": "log_IGRJ16465-4507_18-60.dat",
    "IGR J16479-4514 18-60": "log_IGRJ16479-4514_18-60.dat",
    "IGR J17391-3021 18-60": "log_IGRJ17391-3021_18-60.dat",
    "IGR J17407-2808 18-60": "log_IGRJ17407-2808_18-60.dat",
    "IGR J17544-2619 18-60": "log_IGRJ17544-2619_18-60.dat",
    "IGR J18219-1347 18-60": "log_IGRJ18219-1347_18-60.dat",
    "IGR J18410-0535 18-60": "log_IGRJ18410-0535_18-60.dat",
    "IGR J18450-0435 18-60": "log_IGRJ18450-0435_18-60.dat",
    "IGR J18483-0311 18-60": "log_IGRJ18483-0311_18-60.dat",
    "VelaX-1 18-60": "log_VelaX-1_18-60.dat"
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


def filter_data(data_file):
    return list(filter(
        lambda row: (row["OAA"] < 12) and (row["exposure"] > 500),
        import_data(DATA_FILES[data_file])
    ))


def rebin(rows, bin_size=20):
    for bin_key, bin in itertools.groupby(rows, lambda row: (row["MJD"] // bin_size) * bin_size):
        rows = tuple(bin)
        flux, error = rebin_error(rows)
        yield (bin_key, flux, error)


def rebin_error(rows, keys=("flux", "flux_error")):
    flux = []
    error = []
    for row in rows:
        flux.append(row[keys[0]])
        error.append(row[keys[1]])
    return weighted_av(flux, error), weighted_av_error(flux, error)


def plot_rebin(data):
    x, y, errors = zip(*rebin(data))
    x_shift = empty(len(x))
    for i in range(len(x)):
        x_shift[i] = x[i] - x[0]
    pylab.errorbar(x_shift, y, yerr=errors, fmt=".")
    pylab.xlabel("time, days")
    pylab.ylabel("counts/s")
    pylab.title("VelaX-1 18-60 KeV 20 days binned Lightcurve")
    pylab.show()


def lombscargle(rows, oversample=8., nyquist=10, keys=("MJD", "flux")):
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

    print "period =", period_max

    if plot == True:
        pylab.vlines(period, 0, array(power), color='k', linestyles='solid')
        pylab.xlabel("period, days")
        pylab.ylabel("power")
        pylab.xlim([0, 250])
        pylab.title("VelaX-1 18-60 KeV Periodogram")
        pylab.show()

    return period_max


def fold(rows, period, number_bins=10):
    row_tree = Tree()
    for row in rows:
        row["MJD"] %= period
        lst = row_tree.get(row["MJD"], [])
        lst.append(row)
        row_tree[row["MJD"]] = lst
    #rows.sort(key=itemgetter("MJD"))

    bins = linspace(0, period, number_bins)

    for i in range(len(bins) - 1):
        v = bins[i]
        biggest = bins[i + 1]
        value = itertools.chain.from_iterable(row_tree[v:biggest].values())
        if value:
            flux, error = rebin_error(value)
            yield (v, flux, error)


def plot_fold(data, period_max):
    #x, y, errors = zip(*sorted(fold(data, period_max, 30), key=lambda tup: tup[0]))
    x, y, errors = zip(*fold(data, period_max, 30))
    x = itertools.chain.from_iterable([x, map(lambda x: x + period_max, x)])
    y = itertools.chain.from_iterable([y, y])
    errors = itertools.chain.from_iterable([errors, errors])
    pylab.errorbar(list(x), list(y), yerr=list(errors), fmt=".")
    pylab.xlabel("time, days")
    pylab.ylabel("counts/s")
    pylab.title("VelaX-1 18-60 KeV Folded Lightcurve")
    pylab.show()
    return list(x), list(y), list(errors)


def hratio(high, low):
    x_h, flux_l, err_h = zip(*high)
    x_l, flux_l, err_l = zip(*low)
    ratio = []
    for row in high:
        ratio.append((high[1] - low[1]) / (high[1] + low[1]))
    return x_h, ratio


def hratio_plot(high, low):
    x, y = zip(*hratio(high, low))
    pylab.plot(x, y)
    pylab.xlabel("time")
    pylab.ylabel("hardness ratio")
    pylab.title("VelaX-1 18-60 KeV Folded Lightcurve")
    pylab.show()


def montecarlo(rows):
    for row in rows:
        yield {
            "flux": random.normalvariate(row["flux"], row["flux_error"]),
            "MJD": row["MJD"]
        }


def bootstrap(rows):
    randomrows = []
    for i in len(rows):
        randomrows.append(random.choice(rows))
    for i in list(set(randomrows)):
        yield {
            "flux": row["flux"],
            "MJD": row["MJD"]
        }


def multi_lomb(data, iterations=500):
    p = Pool(4)
    periods = []
    mod_data = []

    for i in range(iterations):
        print i
        mod_data.append(list(montecarlo(data)))

    for item in p.imap(periodogram, mod_data):
        periods.append(item)

    return periods


def histogram(periods, bins=20):
    hist, bins = pylab.histogram(periods, bins)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[: -1] + bins[1:]) / 2
    pylab.bar(center, hist, align='center', width=width)
    pylab.show()


if __name__ == "__main__":
    data = filter_data(data_file="VelaX-1 18-60")
    #data_low = filter_data(data_file="IGR J18027-2016 17-30")
    #data_high = filter_data(data_file="IGR J18027-2016 30-60")
    #data = ["a","b","c","d","e","f"]

    #histogram(multi_lomb(data, iterations=10))

    plot_rebin(data)
    period_max = periodogram(data, plot=True)
    plot_fold(data, period_max)
    #hratio_plot(list(rebin(data_high)), list(rebin(data_low)))
