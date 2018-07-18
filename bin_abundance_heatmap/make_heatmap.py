#!/usr/bin/env python


import sys
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt


# load in abundannce table as pandas dataframe
z=pd.read_csv("bin_abundances.tab", sep='\t', index_col=0)
z.columns.name

# sort by column name
z=z.reindex_axis(sorted(z.columns), axis=1)

# remove all 0 rows
z = z[(z.T != 0).any()]

# standardize columns by total sum in each column
z = z.div(z.sum(axis=0), axis=1)
z=1000000*z

# standardize rows by maximum value in each row
#z = z.div(z.max(axis=1), axis=0)




# remove uninteresting bins
#boi = ["bin.9", "bin.43", "bin.63", "bin.algae"]
#for index in z.index.values:
#        if index not in boi:
#                z=z.drop(index)




# add colored x-labels
lut=[]
for sample in z.columns.values:
        if "2014-09" in sample: lut.append('r')
        elif "2015-06" in sample: lut.append('b')
        elif "2016-02" in sample: lut.append('y')
        elif "2017-02" in sample: lut.append('g')
	else: lut.append('w')


# log standardize:
z=z+0.01
z = np.log(z)

# make heat map
sns.set(font_scale=1)
g = sns.clustermap(z, figsize=(14,8), col_colors=lut, col_cluster=False, cmap="Blues")
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)


plt.show()

