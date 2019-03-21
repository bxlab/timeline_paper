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


def draw_clustermap(filename, ax, c):
	folder = os.path.dirname(os.path.realpath(__file__))
	im = plt.imread(get_sample_data(folder + "/" + filename))
	ax.imshow(im)
	ax.xaxis.set_visible(False)
	ax.set_yticklabels([])
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.yaxis.set_ticks_position('none')



##################   START SCRIPT     ######################


# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)

sns.set_palette("colorblind")
fig = plt.figure(figsize=(10, 6))
colors=["gold", "cyan", "royalblue", "magenta"]
title_font=16

os.system("python make_heatmap.py long_contigs_5kb.tab 9 10 "+" ".join(colors))

ax = fig.add_axes([0.03,  0.08,   0.48,   0.86])
draw_clustermap("heatmap_log.png", ax, colors)
ax.annotate("A", xy=(0.05, 1.02), xycoords="axes fraction", fontsize=title_font)
ax.set_ylabel("Co-assembled contigs (5kb+)", labelpad=-10)
ax.set_title("Contig abundance (log)", fontsize=title_font)


ax = fig.add_axes([0.51,  0.08,   0.48,   0.86])
draw_clustermap("heatmap_std.png", ax, colors)
ax.annotate("B", xy=(0.05, 1.02), xycoords="axes fraction", fontsize=title_font)
ax.set_title("Contig abundance (normalized)", fontsize=title_font)



print "making legend..."
ax = fig.add_axes([0.0, 0.0, 1, 0.05])
ax.axis("off")
legend_elements = [Patch(facecolor=colors[0], edgecolor='k', label='2014', linewidth=1),
        Patch(facecolor=colors[1], edgecolor='k', label='2015', linewidth=1),
        Patch(facecolor=colors[2], edgecolor='k', label='2016', linewidth=1),
        Patch(facecolor=colors[3], edgecolor='k', label='2017', linewidth=1)]
ax.legend(handles=legend_elements, loc="lower center", frameon=True,
	framealpha=1, facecolor='w', ncol=4, columnspacing=1, handlelength=1, 
	prop={'size': 12}, fontsize=12)


#plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)	
plt.savefig("figure_S5.png", dpi=600)



