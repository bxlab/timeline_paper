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

df=pd.read_csv(sys.argv[1], sep='\t', index_col=0)
#df = df.div(df.sum(axis=0), axis=1)
#df*=1000000

df = df.div(df.max(axis=1), axis=0)

for path in df.index.values:
	df.rename({path: " | ".join(path.split(";")[1:])}, axis='index', inplace=True)

# add colored x-labels
lut=[]
for sample in df.columns.values:
        if "2014-09" in sample: lut.append(sys.argv[4])
        elif "2015-06" in sample: lut.append(sys.argv[5])
        elif "2016-02" in sample: lut.append(sys.argv[6])
        elif "2017-02" in sample: lut.append(sys.argv[7])
	else: lut.append('w')

# make heat map
sns.set(font_scale=0.6)
x=float(sys.argv[2])
y=float(sys.argv[3])
g = sns.clustermap(df, figsize=(x,y), col_colors=lut, col_cluster=True, yticklabels=True, xticklabels=False, cmap="magma", vmin=0, vmax=1)
#plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6)
plt.subplots_adjust(left=0, right=1, top=0.99, bottom=0.01)

w_ratio=0.4
h_adjust=0.12
w_adjust=-0.05

hm = g.ax_heatmap.get_position()
xden = g.ax_col_dendrogram.get_position()
yden = g.ax_row_dendrogram.get_position()
col = g.ax_col_colors.get_position()
legend = g.cax.get_position()

g.ax_heatmap.set_position([hm.x0+w_adjust, hm.y0+h_adjust, hm.width*w_ratio, hm.height])
g.ax_col_colors.set_position([col.x0+w_adjust, col.y0+h_adjust, col.width*w_ratio, col.height*0.7])
g.ax_col_dendrogram.set_position([xden.x0+w_adjust, xden.y0-col.height*0.3+h_adjust, xden.width*w_ratio, xden.height*0.5])
g.ax_row_dendrogram.set_position([yden.x0+yden.width*0.3+w_adjust, yden.y0+h_adjust, yden.width*0.7, yden.height])
g.cax.set_position([legend.x0+0.05+w_adjust, legend.y0*0.89+h_adjust, legend.width, legend.height*0.8])

plt.savefig("cluster_pathways.png", dpi=600)
#plt.show()


