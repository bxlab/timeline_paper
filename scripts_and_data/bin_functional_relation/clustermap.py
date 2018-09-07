#!/usr/bin/env python
#./clustermap.py taxonomy.tab bin_funct_annotations/bin.*
import sys
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

taxonomy={}
for line in open(sys.argv[1]):
	if line.startswith("BIN"): continue
	cut=line.strip().split("\t")
	tax=cut[1]
	bin=".".join(cut[0].split(".")[:2])
	taxonomy[bin]=tax.split(";")[2].split("/")[0]


data={}
for f1 in sys.argv[2:]:
	bin=".".join(f1.split("/")[-1].split(".")[:2])
	print bin

	if bin not in data: data[bin]={}
	for line in open(f1):
		cut=line.strip().split("\t")
		info=cut[8]
		if "product=" not in info: continue
		for f in info.split(";"):
			c=f.split("=")
			if c[0]=="product": gene=c[1]
		data[bin][gene]=1



df = pd.DataFrame.from_dict(data)
df=df.fillna(0)


colors={"Actinobacteria": 'b',
"Bacteroidetes": 'r',
"Cyanobacteria": 'y',
"Gammaproteobacteria": 'k',
"Halobacteria": 'g',
"Nanohaloarchaea": 'c',
"unclassified Archaea": 'w',
"Deltaproteobacteria": 'm',
"Chlorophyta": 'w',
"Betaproteobacteria": 'w'}



lut=[]
for bin in df.columns.values:
	tax = taxonomy[bin]
	lut.append(colors[tax])





print "plotting..."
#sns.heatmap(df)
#sns.clustermap(df, col_cluster=True, col_colors=lut, row_colors=lut, robust=True)
sns.clustermap(df, col_colors=lut)

plt.show()





