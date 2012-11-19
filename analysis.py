#!/usr/bin/env python
from __future__ import division
import pylab
from numpy import *
import lomb
import math


MJD, f_flux, f_err, f_expo = loadtxt('log_IGRJ18027-2016_30-60.dat', usecols=(1, 2, 3, 6), unpack=True)
MJD_int = zeros(len(MJD))
MJD_rb = []
f_flux_wav = []


def rebin(MJD):
    for i in range(len(MJD)):
        MJD_int[i] = int(MJD[i])

    startval = MJD_int[0]

    for i in range(len(MJD)):
        MJD_int[i] = MJD_int[i] - startval

    val = 0
    x = []
    delta_x = []
    for i in range(len(MJD)):
        newval = MJD_int[i]
        x.append(f_flux[i])
        delta_x.append(f_err[i])
        if newval != val:
            x.append(f_flux[i])
            delta_x.append(f_err[i])
            f_flux_wav.append(weighted_av(x, delta_x))
            MJD_rb.append(newval)
            x = []
            delta_x = []
            val = newval


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
time = array(MJD_rb)
count = array(f_flux_wav)
fx, fy, nout, jmax, prob = lomb.fasper(time, count, 6., 10)

period_max = 4.57009345794
#period_max = 1 / fx[jmax]
print period_max


def time_f(time):
    for i in range(len(time)):
        time_fold[i] = time[i] % period_max


def count_f():
    for i in range(len(time_fold)):
        for j in range(1, len(time_rb)):
            if time_rb[j] >= time_fold[i]:
                if (time_rb[j - 1] < time_fold[i]):
                    count_rb[j - 1] = count_rb[j - 1] + count[i]

time_fold = zeros(len(time))
time_rb = linspace(0., period_max, num=100)
count_rb = zeros(len(time_rb))
time_f(time)
count_f()


f = open('output_IGRJ18027-2016_30-60.txt', 'w')
for i in range(len(time_rb)):
    print>>f, time_rb[i], count_rb[i]
f.close()

# pylab.plot(time_rb, count_rb, '.')
# pylab.xlabel("time")
# pylab.ylabel("counts")
# pylab.show()
