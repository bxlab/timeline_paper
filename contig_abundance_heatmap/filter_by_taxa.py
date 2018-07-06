#!/usr/bin/env python
import sys


# load taxonomy info
taxa={}
for line in open("contig_taxonomy.tab"):
        cut=line.strip().split("\t")
        #taxa[cut[0]]=cut[1].split(";")
	taxa[cut[0]]=cut[1]


for line in open("contig_abundance.tab"):
	cut=line.strip().split("\t")
	if line[0]=="#": print line.strip()
	if cut[0] not in taxa: continue
	tax=taxa[cut[0]]
	if sys.argv[1] in tax: print line.strip()
