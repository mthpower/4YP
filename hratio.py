#!/usr/bin/env python
import pylab
from numpy import *


time_s, count_s = loadtxt('output_IGRJ18027-2016_17-30.txt', usecols=(0, 1), unpack=True)
time_h, count_h = loadtxt('output_IGRJ18027-2016_30-60.txt', usecols=(0, 1), unpack=True)


def ratio(time_s, time_h, count_s, count_h):
    time = zeros(len(time_h))
    ratio = zeros(len(time_h))
    for i in range(len(time_h)):
        ratio[i] = (count_h[i] - count_s[i]) / (count_h[i] + count_s[i])
    return time, ratio

time, ratio = ratio(time_s, time_h, count_s, count_h)

pylab.plot(time_h, ratio, '.')
pylab.xlabel("time")
pylab.ylabel("hardness ratio")
pylab.title('IGR J18027-2016 Hardness Ratio (30-60 to 17-30 KeV')
pylab.show()
