#!/usr/bin/env python
import sys

for line in open(sys.argv[1]):
	if line[0]=="#": print line.strip()
	else:
		cut=line.strip().split("\t")
		l=int(cut[0].split("_")[3])
		for i in range(1, len(cut)):
			cut[i]=str(float(cut[i])*l)
		print "\t".join(cut)

