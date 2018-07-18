#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
sns.set_palette("colorblind")
plt.rc('font', family='arial')

z=pd.read_csv(sys.argv[1], sep='\t', index_col=0)

# add colored x-labels
lut=[]
for sample in z.columns.values:
        if "2014-09" in sample: lut.append(sys.argv[3])
        elif "2015-06" in sample: lut.append(sys.argv[4])
        elif "2016-02" in sample: lut.append(sys.argv[5])
        elif "2017-02" in sample: lut.append(sys.argv[6])
	else: lut.append('w')

# make heat map
sns.set(font_scale=0.6)
size=float(sys.argv[2])
g = sns.clustermap(z, figsize=(size,size), col_colors=lut, row_colors=lut, col_cluster=True, xticklabels=False, yticklabels=False, cmap="gray")


plt.subplots_adjust(left=0, right=1, top=0.99, bottom=0.01)

ratio=0.6
h_adjust=0.08
w_adjust=-0.05

hm = g.ax_heatmap.get_position()
xden = g.ax_col_dendrogram.get_position()
yden = g.ax_row_dendrogram.get_position()
col = g.ax_col_colors.get_position()
row = g.ax_row_colors.get_position()
legend = g.cax.get_position()

g.ax_heatmap.set_position([hm.x0+w_adjust, hm.y0+h_adjust, hm.width, hm.height])
g.ax_col_colors.set_position([col.x0+w_adjust, col.y0+h_adjust, col.width, col.height])
g.ax_row_colors.set_position([row.x0+w_adjust, row.y0+h_adjust, row.width, row.height])
g.ax_col_dendrogram.set_position([xden.x0+w_adjust, xden.y0+h_adjust, xden.width, xden.height*ratio])
g.ax_row_dendrogram.set_position([yden.x0+w_adjust+yden.width*(1-ratio), yden.y0+h_adjust, yden.width*ratio, yden.height])
g.cax.set_position([legend.x0+0.08+w_adjust, legend.y0*0.93+h_adjust, legend.width, legend.height*0.8])


plt.savefig("uw-unifrac.png", dpi=600)
#plt.show()


