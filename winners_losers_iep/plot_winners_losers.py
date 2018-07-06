#!/usr/bin/env python
import sys, os
from scipy import stats
import numpy as np


#import sample read ct
reads={}
for line in open("sample_read_count.tab"):
	if line[0]=="#": continue
	cut=line.strip().split("\t")
	reads[cut[0]]=int(cut[1])

# import bin taxonomy
taxa={}
for line in open("bin_taxonomy.tab"):
	cut=line.strip().split("\t")
	taxa[cut[0]]=cut[1]

# find winners and losers
verdict={}
fold={}
for line in open("abundance_table.tab"):
	cut=line.strip().split("\t")
	if line.startswith("Genomic"):
		head=cut
	else:
		bin=cut[0]
		pre=[]
		post=[]
		for i in range(1, len(cut)):
			if "2014" in head[i] or "2015" in head[i]:
				pre.append(1000000.0*float(cut[i])/reads[head[i]])
			if "2016" in head[i] or "2017" in head[i]:
				post.append(1000000.0*float(cut[i])/reads[head[i]])

		test = stats.ttest_ind(pre, post)
		fold[bin] = np.mean(post) / np.mean(pre)
		if test[1]>0.05: 
			verdict[bin]="stoic"
			continue
		#if "Halobacteria" not in taxa[bin]: continue

		if np.mean(pre) > np.mean(post):
			verdict[bin]="loser"
		else:
			verdict[bin]="winner"




# load in IEP info
for f in os.listdir("protein_stats"):
	bin='.'.join(f.split('.')[:2])
	if bin not in verdict: continue
	ieps=[]
	for line in open("protein_stats/"+f):
		cut=line.strip().split("\t")
		if len(cut)<5: continue
		if cut[0]=="Pro_id": continue
		iep=float(cut[4])
		ieps.append(iep)
	#print verdict[bin] + "\t" + bin +"\t" + str(np.mean(ieps))
	print bin + "\t" + verdict[bin] + "\t" + str(fold[bin])










		
