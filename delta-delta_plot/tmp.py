#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt

print "loading 2015-2016"
rain=[]
for line in open("2015_2016.tab"):
	rain.append(float(line.split("\t")[0]))

print "loadong 2016-2017"
after=[]
for line in open("2016_2017.tab"):
        after.append(float(line.split("\t")[0]))

print "plotting"
plt.plot(rain, after, 'ko', alpha=0.1, ms=5)
plt.show()
