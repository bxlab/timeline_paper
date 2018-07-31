#!/usr/bin/env python
print "importing python libraries..."
import sys
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


print "loading in abundannce table as pandas dataframe..."
z=pd.read_csv(sys.argv[1], sep='\t', index_col=0)
z.columns.name

# remove all 0 rows
z = z[(z.T != 0).any()]

# standardize columns by total sum in each column
z = z.div(z.sum(axis=0), axis=1)
z=1000000*z


# add colored x-labels
lut=[]
for sample in z.columns.values:
        if "2014" in sample: lut.append(sys.argv[4])
        elif "2015" in sample: lut.append(sys.argv[5])
        elif "2016" in sample: lut.append(sys.argv[6])
        elif "2017" in sample: lut.append(sys.argv[7])
	else: lut.append('w')


print "standardizing by log scale..."
log_df = z+0.01
log_df = np.log(log_df)

print "standardizing columns to the max in each row..."
std_df = z.div(z.max(axis=1), axis=0)


x=float(sys.argv[2])
y=float(sys.argv[3])

print "making log heat map..."
sns.set(font_scale=1)
g = sns.clustermap(log_df, figsize=(x,y), cmap="Blues",
	yticklabels=False, xticklabels=False, col_colors=lut)
plt.savefig("heatmap_log.png", bbox_inches="tight", dpi=600)
plt.clf()


print "making standard heat map..."
g = sns.clustermap(std_df, figsize=(x,y), cmap="Blues", vmin=0, vmax=1,
	yticklabels=False, xticklabels=False, col_colors=lut)
plt.savefig("heatmap_std.png", bbox_inches="tight", dpi=600)
plt.clf()





