#!/usr/bin/env python
import sys
from scipy.stats import pearsonr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

fold_variances={}
rearangements={}

for line in open("data.tab"):
	if line[0]=="#": continue
	cut=line.strip().split("\t")
	
	sampleA = cut[0]
	sampleB=cut[1]
	if sampleA not in fold_variances:
		fold_variances[sampleA]={}
		rearangements[sampleA]={}
	
	fold_var = float(cut[2])
	rearrangement = float(cut[3])
	
	fold_variances[sampleA][sampleB]=fold_var
	rearangements[sampleA][sampleB]=rearrangement

df = pd.DataFrame.from_dict(rearangements)
#df = pd.DataFrame.from_dict(fold_variances)

lut=[]
for sample in df.columns.values:
        if "2014" in sample: lut.append('r')
        elif "2015" in sample: lut.append('b')
        elif "2016" in sample: lut.append('g')
        elif "2017" in sample: lut.append('c')
	else: lut.append('w')



#sns.heatmap(df)
sns.clustermap(df, col_cluster=True, col_colors=lut, row_colors=lut, robust=True, cmap="Blues")


plt.show()





