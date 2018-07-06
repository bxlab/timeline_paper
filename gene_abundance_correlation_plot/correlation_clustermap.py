#!/usr/bin/env python
import sys
from scipy.stats import pearsonr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

corr={}
for f1 in sys.argv[1:]:
	sampleA=f1.split(".")[0]
	for f2 in sys.argv[1:]:
		sampleB=f2.split(".")[0]
		
		if sampleA == sampleB:
			if sampleA not in corr: corr[sampleA]={}
			corr[sampleA][sampleB]=1
			continue

		a=[]
		b=[]
		for line in open(f1): a.append(float(line.strip()))
		for line in open(f2): b.append(float(line.strip()))
	
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



#sns.heatmap(df)
sns.clustermap(df, col_cluster=True, col_colors=lut, row_colors=lut, robust=True)


plt.show()





