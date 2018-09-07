#!/usr/bin/env python
print "loading python packages..."
import sys, getopt, os
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.cbook import get_sample_data
import seaborn as sns
sns.set_color_codes()
import operator as op
import numpy as np
import math
from PIL import Image
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import LogNorm


def draw_clustermap(filename, ax, c, taxa):
	folder = os.path.dirname(os.path.realpath(__file__))
	im = plt.imread(get_sample_data(folder + "/" + filename))
	ax.imshow(im)
	ax.xaxis.set_visible(False)
	ax.set_yticklabels([])
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.yaxis.set_ticks_position('none')
	if taxa=="All":
		ax.set_ylabel("Differentially abundant KEGG pathways", labelpad=-10)
		taxa="All taxa"
	ax.set_title(taxa, fontsize=title_font)



##################   START SCRIPT     ######################


# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)

sns.set_palette("colorblind")
fig = plt.figure(figsize=(10, 6))
colors=["gold", "cyan", "royalblue", "magenta"]
title_font=16


print "making functions clustermaps"
os.system("python make_heatmaps.py 5 5.5 "+" ".join(colors))

coordinates = [ [0.03,	0.06,	0.32,	0.88], 
		[0.35,	0.06,	0.32,	0.88], 
		[0.67,	0.06,	0.32,	0.88] ]
labels = ["A", "B", "C"]
for i,taxa in enumerate(["All", "Archaea", "Bacteria"]):
	ax = fig.add_axes(coordinates[i])
	draw_clustermap(taxa+".png", ax, colors, taxa)
	ax.annotate(labels[i], xy=(0.05, 1.01), xycoords="axes fraction", fontsize=16)


print "making legend..."
ax = fig.add_axes([0.0, 0.0, 1, 0.05])
ax.axis("off")
legend_elements = [Patch(facecolor=colors[0], edgecolor='k', label='2014', linewidth=1),
        Patch(facecolor=colors[1], edgecolor='k', label='2015', linewidth=1),
        Patch(facecolor=colors[2], edgecolor='k', label='2016', linewidth=1),
        Patch(facecolor=colors[3], edgecolor='k', label='2017', linewidth=1)]
ax.legend(handles=legend_elements, loc="lower center", framealpha=1, frameon=True, facecolor='w', ncol=4, columnspacing=1, handlelength=1, prop={'size': 12})


#plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)	
plt.savefig("figure_S4.png", dpi=300)
plt.grid()
#plt.show()



