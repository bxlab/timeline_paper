#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

z=pd.read_csv(sys.argv[1], sep='\t', index_col=0)

# add colored x-labels
lut=[]
for sample in z.columns.values:
        if "2014-09" in sample: lut.append('red')
        elif "2015-06" in sample: lut.append('yellow')
        elif "2016-02" in sample: lut.append('cyan')
        elif "2017-02" in sample: lut.append('magenta')
	else: lut.append('w')

# make heat map
sns.set(font_scale=0.6)
size=float(sys.argv[2])
g = sns.clustermap(z, figsize=(size,size), col_colors=lut, col_cluster=True, xticklabels=False, yticklabels=False, cmap="Blues")

plt.savefig("clustermap.png", bbox_inches='tight')
#plt.show()


