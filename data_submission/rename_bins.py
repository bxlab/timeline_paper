#!/usr/bin/env python
import sys, os

for line in open("rename_bins_table.txt"):
	if line[0]=="#": continue
	cut=line.strip().split("\t")
	orig = cut[0]
	gc=cut[1]
	cov=cut[2]
	tax=cut[3].split(";")[-1]

	print orig +"\t"+ "_".join(["SG1",tax,gc,cov])
	cmd = "cp " + orig+ ".fa " + "_".join(["SG1",tax,gc,cov]) + ".fa"
	os.system( cmd )
	
