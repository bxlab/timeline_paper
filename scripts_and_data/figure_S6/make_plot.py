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




def load_bin_data(filename):
	z=pd.read_csv(filename, sep='\t', index_col=0)
	z = z.div(z.sum(axis=0), axis=1)
	z=1000000*z
	z=z+0.01
	df = np.log(z)
	return df


def color_labels(df):
	# add colored x-labels
	colors=["gold", "cyan", "royalblue", "magenta"]
	lut=[]
	for sample in df.columns.values:
	        if "2014" in sample: lut.append(colors[0])
	        elif "2015" in sample: lut.append(colors[1])
	        elif "2016" in sample: lut.append(colors[2])
	        elif "2017" in sample: lut.append(colors[3])
		else: lut.append('w')
	return lut


def load_taxa(filename):
	taxonomy={}
	for line in open(filename):
		cut = line.strip().split("\t")
		bin=cut[1]
		taxa=cut[2]
		taxonomy[bin]=taxa
	return taxonomy	


def filter_by_taxa(data, taxonomy):
	df = pd.DataFrame()
	for bin,row in data.iterrows():
		taxa=taxonomy[bin]
		if "Cyanobacteria" in taxa:
			bin = taxa.split(";")[-1] +" ("+bin+")"
			if "Total" not in df:
				df["Total"] = row
			df[bin] = row
			df["Total"]+=row
	df.sort_index(axis=0, inplace=True)
	out = pd.DataFrame()
	for sample,row in df.iterrows():
		cut = sample.split("-")
		sample = cut[1]+"-"+cut[3]
		out[sample] = row


	out.drop(["Total"], axis=0, inplace=True)
	return out




data = load_bin_data("bin_abundances.tab")
taxonomy = load_taxa("bin_info.txt")
df = filter_by_taxa(data, taxonomy)
lut = color_labels(df)


sns.set(font_scale=1)
g = sns.clustermap(df, figsize=(10,5), col_colors=lut, col_cluster=False, row_cluster=False, yticklabels=["SG1_Halothece_65_4","SG1_Euhalothece_44_8","SG1_Halothece_48_229"], xticklabels=True, vmin=0, cmap="magma")
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

g.ax_heatmap.set_title("Abundance of Cyanobacteria MAGs", y=1.1, fontsize=14)
plt.savefig("figure_S6.png", dpi=600, bbox_inches="tight")


