#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  3body.py
#
#  Copyright 2012 Matthew <matthew@matthew-VirtualBox>
#-----------------------------------------------------------------------------
#imports

import matplotlib
import pylab
import csv
from decimal import Decimal
counts = []
time = []

with open('testdata.csv', 'r') as data:
    for item in csv.DictReader(data, fieldnames=("time", "counts")):
        counts.append(float(item["counts"]))
        time.append(Decimal(item["time"]))

matplotlib.pyplot.scatter(time, counts)
pylab.xlabel("time")
pylab.ylabel("counts")
pylab.show()
