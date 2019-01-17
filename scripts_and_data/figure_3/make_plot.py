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


#standardize each sample to number of reads in each sample (to 10000 reads per sample)
def standardize_otu_table(df, reads):
	df/=df.sum()
	df*=reads
	return df

# group OTUs by taxonomy
def group_otus_by_taxa(df, rank, taxonomy):
	taxonomy_rank=int(rank)
	taxa_otus={}
	for otu in taxonomy.index.values:
	        if len(taxonomy[otu].split(";"))<taxonomy_rank+1:
	                taxa = "Unknown"
	        else:
	                taxa = taxonomy[otu].split(";")[taxonomy_rank].strip().split("_")[-1]
	        if taxa not in taxa_otus: taxa_otus[taxa]=[]
	        taxa_otus[taxa].append(otu)

	for taxa in taxa_otus:
	        for otu in taxa_otus[taxa]:
	                #rename the index in otu table to taxa
	                df=df.rename(index={otu: taxa})
	#collapse rows with same taxonomy
	df=df.groupby(df.index).sum()
	return df


def load_taxa_data(df, taxa):
	timeline={}
        for sample in df.columns.values:
                year=sample.split("-")[2]+'-'+sample.split("-")[3][:2]
                if year not in timeline: timeline[year]=[]
                timeline[year].append(df.at[taxa, sample])
	return timeline


def load_otu_data(inputfile, rank_of_interest):
	df=pd.read_csv(inputfile, sep='\t', header=1, index_col=0)
	taxonomy=df["taxonomy"]
	df = df.drop("taxonomy", 1)
	
	df = standardize_otu_table(df, 100)
	df = group_otus_by_taxa(df, rank_of_interest, taxonomy)
	return df
	
def calculate_pca(df):
	samples=list(df.index.values)

	pca = PCA(n_components=2)
	principalComponents = pca.fit_transform(df)
	principalDf = pd.DataFrame(data = principalComponents, columns = ['component_1', 'component_2'])

	df = df.reset_index()
	finalDf = pd.concat([principalDf, df[["index"]]], axis = 1)
	finalDf["component_1"]/=1000
	finalDf["component_2"]/=1000
	var = pca.explained_variance_ratio_
	return finalDf, var, samples

def contig_correlation_statistics(data):
	comparisons={}
	for s1 in sorted(data):
		for s2 in sorted(data):
			if s1.split("-")[1] == s2.split("-")[1]: continue
			if int(s1.split("-")[1]) > int(s2.split("-")[1]): continue
			values1 = data[s1].tolist()
			values2 = data[s2].tolist()
			test = stats.pearsonr(values1, values2)
			comparison = s1.split("-")[1]+"-"+s2.split("-")[1]
			if comparison not in comparisons: comparisons[comparison]=[]
			comparisons[comparison].append(test[0])
	print "\t".join(["comparison 1", "comparison 2", "t-statistic", "p-value"])
	for i in sorted(comparisons):
		for j in sorted(comparisons):
			if i==j: continue
			if int(i.split("-")[0])>int(j.split("-")[0]): continue
			test = stats.ttest_ind(comparisons[i], comparisons[j])
			print "\t".join([i,j,str(test[0]), str(test[1])])


def convert_wideform_to_longform(dictionary):
	df={"x":[], "y":[]}
	for key in sorted(dictionary):
		for element in dictionary[key]:
			df["x"].append(key.split("-")[0])
			df["y"].append(element)
	return df

def load_and_standardize_rearrangements(filename):
	values=[]; weights=[]
	for line in open(filename):
		if line[0]=="#": continue
		cut=line.strip().split("\t")
		weights.append(float(cut[0]))
		values.append(float(cut[1]))
	maximum=max(weights)
	out = []
	for i, val in enumerate(values):
		weight = int(100*weights[i]/maximum)
		for j in range(weight):
			out.append(val)
	return out, values

def draw_matrix_clustermap(filename, ax, c):
	folder = os.path.dirname(os.path.realpath(__file__))
	im = plt.imread(get_sample_data(folder + "/" + filename))
	ax.imshow(im)
	ax.xaxis.set_visible(False)
	ax.set_yticklabels([])
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.yaxis.set_ticks_position('none') 
	ax.set_title("Unweighted dissimilarity clustering", fontsize=title_font)
	ax.yaxis.labelpad = 0
	ax.set_ylabel("16S rDNA samples", labelpad=-10)


def draw_bin_clustermap(filename, ax, c):
	folder = os.path.dirname(os.path.realpath(__file__))
	im = plt.imread(get_sample_data(folder + "/" + filename))
	ax.imshow(im)
	ax.xaxis.set_visible(False)
	ax.set_yticklabels([])
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.yaxis.set_ticks_position('none')
	ax.set_title("MAG abundance clustering", fontsize=title_font)
	ax.set_ylabel("Metagenome-assembled genomes", labelpad=-0)


def draw_contig_pca(finalDf, var, samples, c):
	ax.set_title("PCA of contig abundance", fontsize=title_font)
	ax.set_xlabel('PC1 (variance explained = '+ str(var[0]*100)[:4]+"%")
	ax.set_ylabel('PC2 (variance explained = '+ str(var[1]*100)[:4]+"%")
	ax.grid()

	colors={}; i=0
	for sample in sorted(samples):
		label = "-".join(sample.split("-")[1:2])
		if label not in colors:
			colors[label]=c[i]
			i+=1
	for i in finalDf.index:
		x = finalDf["component_1"][i]
		y = finalDf["component_2"][i]
		sample = "-".join(finalDf["index"][i].split("-")[1:2])
		color=colors[sample]
		ax.scatter(x, y, c=color, s=50, edgecolors='k', marker="o", linewidth=0.5)
	for spine in ax.spines.values(): spine.set_alpha(0.2)
	ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)



def draw_violins(data, ax, c):
	c.append(c[0])
	samples=["2014", "2015", "2016", "2017"]
	for s, sample in enumerate(samples):
		df={"x":[], "y":[], "z":[]}
		for i in range(len(data["x"])):
			if data["x"][i]==sample:
				df["z"].append("before")
				df["z"].append("after")
				for j in range(2):
					df["x"].append(data["x"][i])
					df["y"].append(data["y"][i])
		for excluded in samples:
			if excluded==sample: continue
			for j in range(2):
				df["x"].append(excluded)
				df["y"].append(-1)
			df["z"].append("before")
			df["z"].append("after")
		sns.violinplot(x="x", y='y', data=df, width=1, linewidth=1.5, ax=ax, split=True, hue="z", palette=["white", "white"])
		#sns.violinplot(x="x", y='y', data=df, width=1, linewidth=1, ax=ax, split=True, hue="z", palette=[c[s], c[s+1]])

	ax.set_xticklabels(["2014>2015", "2015>2016", "2016>2017", "2017>2014"])
	for tick in ax.get_xticklabels(): tick.set_rotation(20)
	ax.set_ylabel("Average rearrangement index")
	ax.set_ylim(0, 1.05)
	ax.set_xlim(-0.7, 3.7)
	ax.legend_.remove()
	ax.set_title("Taxonomic rearrangement", fontsize=title_font)
	for spine in ax.spines.values(): spine.set_alpha(0.2)
	ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
	plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
	plt.setp(ax.lines, color='k')




def draw_signifficance_bars(df, ax):
	h=96
	for s1 in sorted(data):
		for s2 in sorted(data, reverse=True):
			if s1==s2: continue
			if int(s1.split("-")[0])>int(s2.split("-")[0]): continue
			test=stats.ks_2samp(data[s1], data[s2])
			if test.pvalue > 0.01: continue
			elif test.pvalue > 0.001: m='*'
			elif test.pvalue > 0.0001: m='**'
			else: m='***'

			print s1, s2, m, test
			x_st = int(s1.split("-")[0])-2014
			x_fi = int(s2.split("-")[0])-2014
			ax.hlines(y=h, xmin=x_st, xmax=x_fi, linewidth=1, color='k')
			ax.vlines(x=x_st, ymin=h-1, ymax=h, linewidth=1, color='k')
			ax.vlines(x=x_fi, ymin=h-1, ymax=h, linewidth=1, color='k')
			ax.text(0.95*(x_fi+x_st)/2.0, h-0.5, m)
			h+=2.5

	ax.text(-0.2, 45, "P-value:\n*   0.01\n**  0.001\n*** 0.0001", fontsize=12)



##################   START SCRIPT     ######################


# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)

plt.rc('font', family='arial')
sns.set_palette("colorblind")
#sns.set_style("dark")
fig = plt.figure(figsize=(10, 10))
colors=["gold", "cyan", "royalblue", "magenta"]
title_font=16


print "plotting OTU dissimilarity matrix..."
# xmin, ymin, dx, dy
ax = fig.add_axes([0.045, 0.51, 0.44, 0.44])
os.system("python cluster_matrix.py unweighted_unifrac_matrix.tab 4 "+" ".join(colors))
draw_matrix_clustermap("cluster_matrix.png", ax, colors)
ax.annotate("A", xy=(-0.08, 1.01), xycoords="axes fraction", fontsize=20)


print "making MAG clustermap..."
ax = fig.add_axes([0.52, 0.49, 0.48, 0.46])
os.system("python cluster_mags.py bin_abundances.tab 3.9 4 "+" ".join(colors))
draw_bin_clustermap("cluster_mags.png", ax, colors)
ax.annotate("B", xy=(-0.1, 1.01), xycoords="axes fraction", fontsize=20)


print "making contig abundance PCA..."
ax = fig.add_axes([0.08, 0.12, 0.38, 0.38])
df = pd.read_csv("contig_abundances.tab", delimiter="\t", index_col="#contig").T
df=df.div(df.sum(axis=1)*0.000001, axis=0)
finalDf, var, samples = calculate_pca(df)
contig_correlation_statistics(df.T)
draw_contig_pca(finalDf, var, samples, colors)
ax.annotate("C", xy=(-0.18, 1.01), xycoords="axes fraction", fontsize=20)


print "making rearrangment violin plots"
ax = fig.add_axes([0.56, 0.12, 0.4, 0.38])
os.system("python calculate_rearrangements.py 0")
comparisons=["2014-2015", "2015-2016", "2016-2017", "2017-2014"]
data={}; raw={}
for comparison in comparisons:
	data[comparison], raw[comparison] = load_and_standardize_rearrangements(comparison+".tab")
df = convert_wideform_to_longform(data)
draw_violins(df, ax, colors)
draw_signifficance_bars(raw, ax)
ax.annotate("D", xy=(-0.18, 1.01), xycoords="axes fraction", fontsize=20)


print "making legend..."
ax = fig.add_axes([0.0, 0.0, 1, 0.05])
ax.axis("off")
legend_elements = [Patch(facecolor=colors[0], edgecolor='k', label='2014', linewidth=1),
        Patch(facecolor=colors[1], edgecolor='k', label='2015', linewidth=1),
        Patch(facecolor=colors[2], edgecolor='k', label='2016', linewidth=1),
        Patch(facecolor=colors[3], edgecolor='k', label='2017', linewidth=1)]
ax.legend(handles=legend_elements, loc="lower center", framealpha=1, frameon=True, facecolor='w', ncol=4, columnspacing=1, handlelength=1, prop={'size': 16})


#plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)	
plt.savefig("figure_2.png", dpi=600)
plt.savefig("figure_2.eps", dpi=600)
plt.grid()
#plt.show()



