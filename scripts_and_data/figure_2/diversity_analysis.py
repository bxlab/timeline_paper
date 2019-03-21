#!/usr/bin/env python
import sys
import pandas as pd

depths=pd.read_csv("contig_depth.tab", sep='\t', index_col=0)


genes={}
contigs={}
for sample in depths:
	genes[sample]={}
	contigs[sample]={}

ct=1
for line in open("img_annotation.master"):
	if ct%10000==0:
		print ct
	if line[0]=="#":
		continue
	cut=line.strip().split("\t")
	contig=cut[0].split("-")[0]
	if contig not in depths.index:
		continue
	gene = cut[2]

	for sample in depths:
		if contig not in contigs[sample]:
			contigs[sample][contig]=depths.at[contig, sample]
		if gene not in genes[sample]:
			genes[sample][gene]=0
		genes[sample][gene]+=depths.at[contig, sample]
	ct+=1

for sample in depths:
	line=sample[:]
	for threshold in [1,2,4,8,16,32]:
		ct=0
		for gene in genes[sample]:
			if genes[sample][gene]>threshold:
				ct+=1
		#for contig in contigs[sample]:
		#	if contigs[sample][contig]>threshold:
		#		ct+=int(contig.split("_")[3])
		line+="\t"+str(ct)
	print line





