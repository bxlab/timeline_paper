#!/usr/bin/env python


import sys
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# load in abundannce table as pandas dataframe
z=pd.read_csv(sys.argv[1], sep='\t', index_col=0)
z.columns.name


# add colored x-labels
lut=[]
for sample in z.columns.values:
	if "2014" in sample: lut.append('r')
        elif "2015" in sample: lut.append('b')
        elif "2016" in sample: lut.append('y')
        elif "2017" in sample: lut.append('g')
	else: lut.append('w')


# log standardize:
z=z+1
z = np.log(z)


# make heat map
sns.set(font_scale=1)
g = sns.clustermap(z, figsize=(8,8), col_colors=lut, row_colors=lut, col_cluster=False, row_cluster=False)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)



#plt.savefig(sys.argv[3], bbox_inches='tight')
plt.show()


