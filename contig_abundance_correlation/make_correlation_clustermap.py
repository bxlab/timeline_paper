#!/usr/bin/env python
import sys
from scipy.stats import pearsonr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv(sys.argv[1], delimiter="\t", index_col="#contig")

lists={}
for sample in df:
	lists[sample]=df[sample].tolist()

corr = {}
for sampleA in df:
	for sampleB in df:
		if sampleA == sampleB:
			if sampleA not in corr: corr[sampleA]={}
			corr[sampleA][sampleB]=1
			continue
		a=lists[sampleA]
		b=lists[sampleB]
		r,p = pearsonr(a,b)
		if sampleA not in corr: corr[sampleA]={}
		corr[sampleA][sampleB]=r

df = pd.DataFrame.from_dict(corr)

lut=[]
for sample in df.columns.values:
        if "2014" in sample: lut.append('r')
        elif "2015" in sample: lut.append('b')
        elif "2016" in sample: lut.append('g')
        elif "2017" in sample: lut.append('c')
	else: lut.append('w')


sns.clustermap(df, col_cluster=True, col_colors=lut, row_colors=lut)


plt.show()

