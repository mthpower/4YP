#!/usr/bin/env python
from __future__ import division
import itertools
import pylab
from numpy import *
import lomb
import math


#MJD, f_flux, f_err, f_expo = loadtxt('log_GX13+1_30-60.dat', usecols=(1, 2, 3, 6), unpack=True)
#MJD_int = zeros(len(MJD))
MJD_rb = []
f_flux_wav = []
f_err_rb = []


def import_data(file_name):
    """
    Imports text file of format specified by T.Bird, outputs into arrays:
    MDJ, flux, flux_error, OOA, exposure.
    example MJD, flux, flux_error, OOA, exposure = import_data(file_name)
    """
    return loadtxt(file_name, usecols=(1, 2, 3, 5, 6), dtype={
        "names": ("MJD", "flux", "flux_error", "OOA", "exposure"),
        "format": ("f", "f", "f", "f", "f")
    })

bin_size = 5
itertools.regroup(import_data("log_GX13+1_30-60.dat"), lambda row: (row["MJD"] // bin_size) * bin_size)

def rebin(rows):
    """
    Bins input arrays into smaller bins specified by bin_size.
    Input arrays are MJD and flux and flux_error
    Output arrays are rebinned_MJD, rebinned_flux, rebinned_flux_error
    """
    # for i in range(len(MJD)):
    #     MJD_int[i] = int(MJD[i])

    rebinned_MJD = arange(
        MJD[0],
        (MJD[-1] + (bin_size / 100)),
        bin_size
    )





    for i in range(len(rebinned_MJD)):
        rebinned_MJD[i] = rebinned_MJD[i] - MJD[0]

    bin_check = 0
    fluxes = []
    errors = []
    for i in range(len(rebinned_MJD) - 1):
        for j in range(len(MJD)):
            if MJD[j] >= rebinned_MJD[i] and MDJ[j] < rebinned_MJD[i + 1]:
                fluxes.append(flux[j])
                errors.append(flux_error[j])

    yield (rebinned_MJD,foo_bin)

    for MJD_bin, foo_bin in rebin(*args):
        dostuff(MJD_bin)
        dostuff(foo_bin)


    delta_x = []
    for i in range(len(rebinnedMJD))l:
        newval = rebinnedMJD[i]
        x.append(f_flux[i])
        delta_x.append(f_err[i])
        if newval != val:
            x.append(f_flux[i])
            delta_x.append(f_err[i])
            f_err_rb.append(error(delta_x))
            f_flux_wav.append(weighted_av(x, delta_x))
            MJD_rb.append(newval)
            x = []
            delta_x = []
            val = newval


def error(errors):
    sqerrs = zeros(len(errors))
    for i in range(len(errors)):
        sqerrs[i] = errors[i]
    return math.sqrt(math.fsum(sqerrs))


def weighted_av(x, delta_x):
    w = zeros(len(delta_x))
    for i in range(len(delta_x)):
        w[i] = 1 / ((delta_x[i]) ** 2)
    sum_w = math.fsum(w)
    xw = zeros(len(x))
    for i in range(len(x)):
        xw[i] = x[i] * w[i]
    sum_xw = math.fsum(xw)
    wav = (sum_xw / sum_w)
    return wav


rebin(MJD)
time = array(MJD)
count = array(f_flux)
fx, fy, nout, jmax, prob = lomb.fasper(MJD, f_flux, 6., 10)

period_max = 1 / fx[jmax]
print period_max

period = zeros(len(fx))
for i in range(len(fx)):
    period[i] = 1 / fx[i]


# pylab.plot(period, fy, '.')
# pylab.xlabel("period, days")
# pylab.ylabel("power")
# pylab.xlim([0,300])
# pylab.title('GX13+1 30-60 KeV Periodogram ')
# pylab.show()


def time_f(time):
    for i in range(len(time)):
        time_fold[i] = time[i] % (2 * period_max)


def count_f():
    for i in range(len(time_fold)):
        for j in range(1, len(time_rb)):
            if time_rb[j] >= time_fold[i] and (time_rb[j - 1] < time_fold[i]):
                count_rb[j - 1] = count_rb[j - 1] + count[i]


def error_f():
    for i in range(len(time_fold)):
        for j in range(1, len(time_rb)):
            if time_rb[j] >= time_fold[i] and (time_rb[j - 1] < time_fold[i]):
                count_rb[j - 1] = count_rb[j - 1] + count[i]


time_fold = zeros(len(time))
time_rb = linspace(0., period_max)
count_rb = zeros(len(time_rb))
error_rb = zeros(len(time_rb))
time_f(time)
count_f()

# f = open('output_IGRJ18027-2016_17-30.txt', 'w')
# for i in range(len(time_rb)):
#     print>>f, time_rb[i], count_rb[i]
# f.close()

pylab.plot(time_rb, count_rb, '.')
#pylab.errorbar(MJD_rb, f_flux_wav, yerr=f_err_rb, fmt=".")
pylab.xlabel("time")
pylab.ylabel("counts")
pylab.title('GX13+1 30-60 KeV Folded Lightcurve')
pylab.show()