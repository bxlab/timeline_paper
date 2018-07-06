#!/usr/bin/env python2
import sys
import matplotlib.pyplot as plt

data=[]
for line in open(sys.argv[1]):
	cut=line.strip().split('\t')
	if line[0]=="#": 
		print line.strip()
		continue
	if int(cut[0].split("_")[3])<5000: break
	print line.strip()

	#for plotting histogram
	#tot=0
	#for f in cut[1:]:
	#	tot+=float(f)
	#if tot<3000000: data.append(tot)
	
	#if tot>1000000: print line.strip()
	

#plt.hist(data, 150, density=True, facecolor='g', alpha=0.75)
#plt.show()
