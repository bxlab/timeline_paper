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
		pal[i]="skyblue"

	sns.boxplot(data=values, width=0.60, linewidth=1, ax=ax, palette=pal)
	sns.swarmplot(data=values, size=4, edgecolor="black", linewidth=1, ax=ax, palette=pal)

	ax.set_xticklabels(keys)
	for tick in ax.get_xticklabels(): tick.set_rotation(45)
	for spine in ax.spines.values(): spine.set_alpha(0.2)





##################   START SCRIPT     ######################

# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
fig = plt.figure(figsize=(10, 6))
colors=["gold", "cyan", "royalblue", "magenta"]
title_font=16

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
	if i>1: i=i+2
	ax = fig.add_subplot(2,4,i+1)
	if i<2: ax.get_xaxis().set_ticks([])
	ax.annotate(letter, xy=(-0.15, 1.02), xycoords="axes fraction", fontsize=title_font)
	if taxa=="Cytophagia": name="Bacteroidetes"
	else: name=taxa
	ax.set_title(name, fontsize=title_font)
	data=taxa_data[taxa]
	draw_boxplots(data, taxa, ax, colors)
	if i<2: ax.get_xaxis().set_ticks([])


print "making domain abundance plot"
taxa_of_interest=["Archaea"]
df = load_otu_data("otu_table.tab", 0)
taxa_data={}
for taxa in taxa_of_interest:
        taxa_data[taxa]=load_taxa_data(df, taxa)

letter="E"
ax = fig.add_subplot(1,2,2)
ax.annotate(letter, xy=(-0.05, 1.02), xycoords="axes fraction", fontsize=title_font)
ax.set_title(taxa, fontsize=title_font)
data=taxa_data[taxa]
draw_boxplots(data, taxa, ax, colors)
if i==0: ax.set_ylabel("Relative taxa abundance (%)")



plt.tight_layout()
plt.savefig("figure_S3.png", dpi=300)
#plt.show()



