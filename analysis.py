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
import random
from multiprocessing import Pool

DATA_FILES = {
    "IGR J18027-2016 17-30": "log_IGRJ18027-2016_17-30.dat",
    "IGR J18027-2016 30-60": "log_IGRJ18027-2016_30-60.dat",
    "GX13+1_30-60": "log_GX13+1_30-60.dat"
}


# OBJECTS = (
#     {
#         "title": "J18027-2016",
#         "ranges": (
#             {
#                 "bands": "30-60",
#                 "file_name": "log_IGRJ18027-2016_17-30.dat"
#             }
#         # ...
#         )
#     }
#     # ...
# )


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
    x_shift = zeros(len(x))
    for i in range(len(x)):
        x_shift[i] = x[i] - x[0]
    pylab.errorbar(x_shift, y, yerr=errors, fmt=".")
    pylab.xlabel("time")
    pylab.ylabel("counts")
    pylab.title("IGRJ18027 2016 17-30 KeV 1 Month Lightcurve")
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

    print "period =", period_max

    if plot == True:
        pylab.vlines(period, 0, array(power), color='k', linestyles='solid')
        pylab.xlabel("period, days")
        pylab.ylabel("power")
        pylab.xlim([0, 40])
        pylab.title("IGRJ18027 2016_17-30 KeV Periodogram")
        pylab.show()

    return period_max


def fold(rows, period, number_bins=10):
    for row in rows:
        row["MJD"] %= period
    rows.sort(key=itemgetter("MJD"))

    MJD_bins = {k: list() for k in linspace(0, period, number_bins)[:-1]}  # dropping the last element

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
    #x, y, errors = zip(*sorted(fold(data, period_max, 30), key=lambda tup: tup[0]))
    x, y, errors = zip(*fold(data, period_max, 30))
    x = itertools.chain.from_iterable([x, map(lambda x: x + period_max, x)])
    y = itertools.chain.from_iterable([y, y])
    errors = itertools.chain.from_iterable([errors, errors])
    pylab.errorbar(list(x), list(y), yerr=list(errors), fmt=".")
    pylab.xlabel("time")
    pylab.ylabel("counts")
    pylab.title("IGRJ18027 2016 17-30 KeV Folded Lightcurve")
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
    pylab.title("IGRJ18027 2016 17-30 KeV Folded Lightcurve")
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


# def lombscargle_multi(data):
#     period = periodogram(montecarlo(data))
#     return period


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
    data = filter_data(data_file="IGR J18027-2016 17-30")
    #data_low = filter_data(data_file="IGR J18027-2016 17-30")
    #data_high = filter_data(data_file="IGR J18027-2016 30-60")
    #data = ["a","b","c","d","e","f"]

    histogram(multi_lomb(data))

    #plot_rebin(data)
    #period_max = periodogram(data)
    #plot_fold(data, period_max)
    #hratio_plot(list(rebin(data_high)), list(rebin(data_low)))
