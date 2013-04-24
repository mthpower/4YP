#!/usr/bin/env python
from __future__ import division
import itertools
import pylab
import numpy as np
import scipy.signal
import lomb
import math
from operator import mul
from bintrees import FastBinaryTree as Tree
from six.moves import filter, zip
import random
from multiprocessing import Pool

from nose.tools import eq_


def quadrature(dx):
    sqerrs = map(lambda dx: dx ** 2, dx)
    return math.sqrt(math.fsum(sqerrs))


def fractional(x, dx):
    sqfracerrs = map(lambda x, dx: (dx / x) ** 2, x, dx)
    mult = reduce(mul, x)
    return mult * math.sqrt(math.fsum(sqfracerrs))


def weighted_av(x, dx):
    weights = map(lambda x: 1 / (x ** 2), dx)
    xw = map(mul, x, weights)
    inverted_dx = map(lambda x: 1 / x, dx)
    dU = quadrature(inverted_dx)
    return math.fsum(xw) / math.fsum(weights), dU / math.fsum(weights)


# def weighted_av_error(x, dx):
#     def fractxy(x, dx, y):
#         return y * (math.sqrt((dx / x) ** 2))
#     w = map(lambda x: 1 / (x ** 2), dx)
#     q = map(mul, x, w)
#     dq = map(fractxy, x, dx, q)
#     return fractxy(math.fsum(q), quadrature(dq), weighted_av(x, dx))


def import_data(file_name):
    """
    Imports text file of format specified by T.Bird, outputs into arrays:
    MJD, flux, flux_error, OOA, exposure.
    example MJD, flux, flux_error, OOA, exposure = import_data(file_name)
    """
    return np.loadtxt(file_name, usecols=(1, 2, 3, 5, 6), dtype={
        "names": ("MJD", "flux", "flux_error", "OAA", "exposure"),
        "formats": ("f", "f", "f", "f", "f")
    })

def import_data_asm(file_name):
    """
    Imports text file of format specified by T.Bird, outputs into arrays:
    MJD, flux, flux_error, OOA, exposure.
    example MJD, flux, flux_error, OOA, exposure = import_data(file_name)
    """
    return np.loadtxt(file_name, usecols=(0, 1, 2), dtype={
        "names": ("MJD", "flux", "flux_error"),
        "formats": ("f", "f", "f")
    })


def filter_data(file_name):
    return list(filter(
        lambda row: (row["OAA"] < 12) and (row["exposure"] > 500),
        import_data(file_name)
    ))


def rebin(rows, bin_size=20):
    """Rows must be sorted by MJD"""
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
    return weighted_av(flux, error)  # , weighted_av_error(flux, error)


def plot_rebin(data):
    x, y, errors = zip(*rebin(data))
    x_shift = np.empty(len(x))
    for i in range(len(x)):
        x_shift[i] = x[i] - x[0]
    #pylab.errorbar(x_shift, y, yerr=errors)
    pylab.errorbar(x_shift, y, yerr=errors, fmt=".")
    pylab.xlabel("time, days")
    pylab.ylabel("counts/s")
    pylab.title("CenX-3 18-60 KeV 20 days binned Lightcurve")
    pylab.show()


def lombscargle(rows, oversample=6, nyquist=5, keys=("MJD", "flux")):
    x = []
    y = []
    for row in rows:
        x.append(row[keys[0]])
        y.append(row[keys[1]])

    fx, fy, nout, jmax, prob = lomb.fasper(np.array(x), np.array(y), oversample, nyquist)

    return (zip(fx, fy), jmax)


def lombscargle_scipy(rows, frequencies, keys=("MJD", "flux")):
    x = np.empty(len(rows), dtype=np.float64)
    y = np.empty(len(rows), dtype=np.float64)
    for i, row in enumerate(rows):
        x[i] = row[keys[0]]
        y[i] = row[keys[1]]

    pgram = scipy.signal.lombscargle(x, y, frequencies)
    return frequencies, pgram


def periodogram(data, plot=False):
    data, jmax = lombscargle(data)

    # calculate the reciprocal of the x dimension
    data = list(itertools.starmap(lambda x, y: (1 / x, y), data))
    period, power = zip(*data)
    period_max = data[jmax][0]

    print "period =", period_max, "days"

    if plot:
        pylab.vlines(period, 0, np.array(power), color='k', linestyles='solid')
        pylab.xlabel("period, days")
        pylab.ylabel("power")
        pylab.xlim([2, 200])
        pylab.title("CenX-3 18-60 KeV Periodogram")
        pylab.show()

    return period_max


def periodogram_scipy(data, frequencies, plot=False):
    freqs, pgram = lombscargle_scipy(data, frequencies)
    periods = map(lambda x: (2 * math.pi) / x, freqs)

    if plot:
        pylab.vlines(periods, 0, np.array(pgram), color='k', linestyles='solid')
        pylab.xlabel("period, days")
        pylab.ylabel("power")
        pylab.title("CenX-3 18-60 KeV Periodogram")
        pylab.show()


def fold(rows, period, number_bins=20, flatten=True):
    row_tree = Tree()
    for row in rows:
        mod_mjd = row["MJD"] % period
        lst = row_tree.get(mod_mjd, [])
        lst.append(row)
        row_tree[mod_mjd] = lst

    bins = np.linspace(0, period, number_bins + 1)
    for i in range(len(bins) - 1):
        v = bins[i]
        biggest = bins[i + 1]
        value = list(itertools.chain.from_iterable(row_tree[v:biggest].values()))
        if len(value):
            if not flatten:
                yield (v, bin_fluxes(value))
            else:
                flux, error = rebin_error(value)
                yield (v, flux, error)


def plot_fold(data, period_max):
    #x, y, errors = zip(*sorted(fold(data, period_max, 30), key=lambda tup: tup[0]))
    x, y, errors = zip(*fold(data, period_max, 30))
    x = itertools.chain.from_iterable([x, map(lambda x: x + period_max, x)])
    y = itertools.chain.from_iterable([y, y])
    errors = itertools.chain.from_iterable([errors, errors])
    x = map(lambda x: x / period_max, x)
    pylab.errorbar(list(x), list(y), yerr=list(errors))
    #pylab.errorbar(list(x), list(y), yerr=list(errors), fmt=".")
    pylab.xlabel("phase, cycles")
    pylab.ylabel("counts/s")
    pylab.title("CenX-3 18-60 KeV Folded Lightcurve")
    pylab.show()
    return list(x), list(y), list(errors)


def _inner_pdm(rows, period):
    #folded = list(fold(rows, period, flatten=False))
    bin, fluxes = zip(*fold(rows, period, 30, flatten=False))
    M = len(bin)
    variance = []
    n = []
    for flux in fluxes:
        variance.append(np.array(flux).var())
        n.append(len(flux))
    return s_variance(variance, n, M)


def calculate_sigma(rows):
    lc_fluxes = []
    for row in rows:
        lc_fluxes.append(row["flux"])
    return np.array(lc_fluxes).var()


def pdm(rows, frequencies):
    periods = map(lambda x: (2 * math.pi) / x, frequencies)
    sigma = calculate_sigma(rows)
    for period in periods:
        print "folding on" , period 
        yield (_inner_pdm(rows, period)) / sigma


def bin_fluxes(rows):
    fluxes = []
    for row in rows:
        fluxes.append(row["flux"])
    return fluxes


def s_variance(variance, n, M):
    A = map(lambda var, n: (n - 1) * var, variance, n)
    return math.fsum(A) / (math.fsum(n) - M)


def plot_pdm(rows, frequencies):
    periods = map(lambda x: (2 * math.pi) / x, frequencies)
    pylab.plot(periods, list(pdm(rows, frequencies)))
    pylab.xlabel("period, days")
    pylab.ylabel("dispersion")
    pylab.title("IGRJ18027-2016 18-60 KeV Phase Dispersion")
    pylab.show()


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
    data = filter_data(file_name="project_data/log_IGRJ16465-4507_18-60.dat")
    #data_asm = import_data_asm(file_name="ASM/ASMCenX3.txt")
    #data_low = filter_data(data_file="IGR J18027-2016 17-30")l
    #data_high = filter_data(data_file="IGR J18027-2016 30-60")
    #data = ["a","b","c","d","e","f"]

    #histogram(multi_lomb(data, iterations=10))
    #plot_rebin(data_asm)
    plot_rebin(data)
    period_max = periodogram(data, plot=False)
    #period_max = periodogram(data_asm, plot=False)
    plot_fold(data, period_max)
    #plot_fold(data_asm, period_max)
    #print pdm(data, arange(10.1, 10.4, 0.1))
    # data = []
    # for i in np.arange(1, 100, 0.1):
    #     data.append({
    #         "MJD": i,
    #         "flux": math.sin(i / 3),
    #         "error": 1,
    #     })

    #periodogram_scipy(data, np.arange(0.02, 5, 0.0001), plot=True)

    #plot_pdm(data, np.arange(0.18, 4, 0.005))


def test_pdm():
    data = filter_data(file_name="log_VelaX-1_18-60.dat")
    print _inner_pdm(data, 10.1)
    print _inner_pdm(data, 10.1)

    sigma = calculate_sigma(data)

    bin_sigmas = []
    for period in np.arange(10.0, 10.4, 0.1):
        bin_sigmas.append(_inner_pdm(data, period) / sigma)

    print bin_sigmas

    sigma = calculate_sigma(data)
    bin_sigmas = []
    for period in np.arange(10.2, 10.4, 0.1):
        bin_sigmas.append(_inner_pdm(data, period) / sigma)

    print bin_sigmas

    #plot_fold(data, period_max)
    #hratio_plot(list(rebin(data_high)), list(rebin(data_low)))


def test_pdm_sin():
    data = []
    for i in np.arange(1, 100, 0.1):
        data.append({
            "MJD": i,
            "flux": math.sin(i / 3),
            "error": 1,
        })

    eq_(list(pdm(data, np.arange(19, 19.3, 0.1))),
        [0.0086289659844807666, 0.017522407256532074, 0.030554646990128276, 0.047455228134305749]
        )
