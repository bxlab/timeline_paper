#!/usr/bin/env python
import sys
from scipy.stats import pearsonr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

df = pd.read_csv(sys.argv[1], delimiter="\t", index_col="Category")
df = df[~df.index.duplicated(keep='first')]

lut=[]; samples=[]
for sample in df.columns.values:
        if "2014" in sample: 
		lut.append('r')
		samples.append("2014-"+sample.split("-")[3])
        elif "2015" in sample:
		lut.append('g')
		samples.append("2015-"+sample.split("-")[3])
        elif "2016" in sample: 
		lut.append('b')
		samples.append("2016-"+sample.split("-")[3])
        elif "2017" in sample:
		lut.append('c')
		samples.append("2017-"+sample.split("-")[3])
	else: lut.append('w')


# standardize columns by total sum in each column
df = df.div(df.sum(axis=0), axis=1); df=1000000*df


# test significance
tot=0
for gene in df.index.values:
	t2014=[];t2015=[];t2016=[];t2017=[]
	for sample in df.columns.values:
		if "2014" in sample:
			t2014.append(df[sample][gene])
		if "2015" in sample:
			t2015.append(df[sample][gene])
		if "2016" in sample:
			t2016.append(df[sample][gene])
		if "2017" in sample:
			t2017.append(df[sample][gene])
	anova = stats.f_oneway(t2014, t2015, t2016, t2017)
	if anova[1]>0.001:
		keep=False
	else:
		keep=True

	if keep==True:
		fold = np.mean(t2016)/np.mean(t2015)
		print  str(fold)[:4] + '\t' + gene
		tot+=1
	else:
		df=df.drop(gene)

print str(tot) + " signifficant functional groups found!"


# log standardize
#df=df+1
#df = np.log(df)

# standardize rows by maximum value in each row
df = df.div(df.max(axis=1), axis=0)
g = sns.clustermap(df, col_cluster=True, col_colors=lut, cmap="Blues", figsize=(10,8), vmin=0, vmax=1, xticklabels=samples)

plt.subplots_adjust(left=0.05, bottom=0.1, right=0.7, top=0.95)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.show()











