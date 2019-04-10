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
from PIL import Image

def clear_ax(ax):
	ax.xaxis.set_visible(False)
	ax.set_yticklabels([])
	ax.set_yticklabels([])
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.yaxis.set_ticks_position('none')
	


def insert_png(filename, ax):
	folder = os.path.dirname(os.path.realpath(__file__))
	im = plt.imread(get_sample_data(folder + "/" + filename))
	ax.imshow(im)
	ax.set_title(filename.split(".")[0], fontsize=title_font)


def manual_png(img, x, y):
	print "plotting "+img
	folder = os.path.dirname(os.path.realpath(__file__))
	im = Image.open(folder + "/" + img)

	plot_h = fig.bbox.ymax*1
	plot_w = fig.bbox.xmax*0.9

	h = im.size[1]; w = im.size[0]
	resize_ratio = min(plot_h/h, plot_w/w)
	h*=resize_ratio; w*=resize_ratio

	im = im.resize((int(w), int(h)), Image.ANTIALIAS)
	im = np.array(im).astype(np.float) / 255
	plt.figimage(im, x, y)
	

def add_labels(labels):
	ax = fig.add_axes([0, 0, 1, 1])
	ax.text(0.015, 0.93, labels[0], fontsize=26)
	ax.text(0.015, 0.58, labels[1], fontsize=26)
	ax.axis("off")
	

##################   START SCRIPT     ######################


# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
fig = plt.figure(figsize=(8,8), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])

# xmin, ymin, dx, dy
x_max = fig.bbox.xmax
y_max = fig.bbox.ymax
print x_max, y_max
manual_png("subfig_A.png", x_max*0.06, y_max*0.65)
manual_png("subfig_B.png", x_max*0.06, y_max*0.03)


add_labels(["A","B"])

plt.savefig("figure_S1.png", dpi=300)
#plt.show()
