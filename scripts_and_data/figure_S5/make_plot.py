#!/usr/bin/env python
print "loading python packages"
import sys, getopt, os
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
sns.set_color_codes()
import operator as op
import numpy as np
import math
from PIL import Image

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
	

def convert_lists_to_df(x,y,labels):
	df = pd.DataFrame({'labels': labels,'x': x,'y': y})
	return df


def draw_boxplots(data, taxa, ax, col):
	keys=[]; values=[]; pal={}
	for i, key in enumerate(sorted(data)):
		keys.append(key)
		values.append(data[key])
                if key=="2016-02": c=col[0]
                if key=="2016-07": c=col[1]
                if key=="2016-10": c=col[2]
                if key=="2017-02": c=col[3]
                pal[i]=c

	sns.boxplot(data=values, width=0.60, linewidth=1, ax=ax, palette=pal)
	sns.swarmplot(data=values, size=4, edgecolor="black", linewidth=1, ax=ax, palette=pal)

	ax.set_xticklabels(keys)
	for tick in ax.get_xticklabels(): tick.set_rotation(30)
	for spine in ax.spines.values(): spine.set_alpha(0.2)
	ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)


def get_max_min_in_dict(dictionary):
	maxs=[]; mins=[]
	for key in dictionary:
		maxs.append(max(dictionary[key]))
		mins.append(min(dictionary[key]))
	return max(maxs), min(mins)


def draw_signifficance_bars(df, ax, inc_mod=1):
	max_val,min_val = get_max_min_in_dict(df)
	h = (max_val-min_val)*1.1 + min_val
	inc = (h-min_val)*0.05*inc_mod

	for x_st,s1 in enumerate(sorted(df)):
		for x_fi,s2 in enumerate(sorted(df)):
			s1_tot=float(s1.split("-")[0]) + float(s1.split("-")[0])/10
			s2_tot=float(s2.split("-")[0]) + float(s2.split("-")[0])/10
			if s1>=s2: continue
			test=stats.ttest_ind(df[s1], df[s2])
			if test.pvalue > 0.01: continue
			elif test.pvalue > 0.001: m='*'
			elif test.pvalue > 0.0001: m='**'
			else: m='***'
			ax.hlines(y=h, xmin=x_st, xmax=x_fi, linewidth=0.5, color='k')
			ax.text((x_fi+x_st)/2.0, h-inc/1.5, m, ha='center', fontsize=10)
			h+=inc



##################   START SCRIPT     ######################

# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
fig = plt.figure(figsize=(5, 5))
title_font=14
colors=["royalblue", "cyan", "magenta", "red"]


main = fig.add_subplot(111)
#main.set_title("Taxa abundance recovery post-rain", fontsize=title_font, y=1.1)
main.set_ylabel("Relative taxa abundance (%)", labelpad=20)
main.xaxis.set_visible(False)
main.set_yticklabels([])
for spine in main.spines.values(): spine.set_visible(False)
main.yaxis.set_ticks_position('none')

print "making phyla abundance plot"
taxa_of_interest=["Cyanobacteria","Chloroplast","Cytophagia","Halobacteria"]	
df = load_otu_data("otu_table.tab", 2)
taxa_data={}
for taxa in taxa_of_interest: 
	taxa_data[taxa]=load_taxa_data(df, taxa)

labels=["A","B","C","D"]
for i, taxa in enumerate(taxa_data):
	letter=labels[i]
	#if i>1: i=i+2
	ax = fig.add_subplot(2,2,i+1)
	if i<2: ax.get_xaxis().set_ticks([])
	ax.annotate(letter, xy=(-0.15, 1.03), xycoords="axes fraction", fontsize=title_font)
	if taxa=="Cytophagia": name="Bacteroidetes"
	elif taxa=="Chloroplast": name="Chlorophyta"
	else: name=taxa
	ax.set_title(name, fontsize=title_font)
	data=taxa_data[taxa]
	draw_boxplots(data, taxa, ax, colors)
	draw_signifficance_bars(data, ax)
	if i<2: ax.get_xaxis().set_ticklabels(["","","",""])


plt.tight_layout(rect=[-0.04, -0.0, 1, 1])
plt.savefig("figure_S5.png", dpi=600)



