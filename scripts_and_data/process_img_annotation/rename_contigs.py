#!/usr/bin/env python
import sys

names={}
for line in open(sys.argv[1]):
	cut=line.strip().split("\t")
	names[cut[1]]=cut[0]

for line in open(sys.argv[2]):
	cut=line.strip().split("\t")
	contig=cut[0][:17]
	id=cut[0][17:]
	contig = names[contig]
	cut[0]=contig+"-"+id
	
	if int(cut[0].split("_")[3])<1000: continue

	print "\t".join(cut)
