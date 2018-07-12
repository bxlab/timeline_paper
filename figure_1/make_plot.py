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
	

def load_pca_coord(filename):
	pc1=[]; pc2=[]; labels=[]	
	for i,line in enumerate(open(filename)):
		cut=line.strip().split("\t")
		if i==4:
			expl_1=float(cut[0])
			expl_2=float(cut[1])
		if i>8:
			if len(cut)<10: continue
			pc1.append(float(cut[1]))
			pc2.append(float(cut[2]))
			labels.append("-".join(cut[0].split("-")[2:4])[:7])
	return expl_1,expl_2,pc1,pc2,labels


def convert_lists_to_df(x,y,labels):
	df = pd.DataFrame({'labels': labels,'x': x,'y': y})
	return df


def draw_boxplots(data, taxa, ax, position):
	keys=[]; values=[]; pal={}
	for i, key in enumerate(sorted(data)):
		keys.append(key)
		values.append(data[key])
		if key=="2014-09" or key=="2015-06":
			pal[i]='yellow'
		else:
			pal[i]='cyan'

	sns.boxplot(data=values, width=0.60, linewidth=1, ax=ax, palette=pal)
	sns.swarmplot(data=values, size=3, edgecolor="black", linewidth=1, ax=ax, palette=pal)

	ax.set_xticklabels(keys)
	for tick in ax.get_xticklabels(): tick.set_rotation(45)

	for spine in ax.spines.values(): spine.set_alpha(0.2)

	if taxa=="Cytophagia": taxa="Bacteroidetes"
	if position==1 or position==0: ax.set_title(taxa)


def draw_pcoa(PC1_expl, PC2_expl, df, ax):
	plots={}
	for i,label in enumerate(df["labels"]):
		if "2014" in label: c='red'
		if "2015" in label: c='yellow'
		if "2016" in label: c='cyan'
		if "2017" in label: c='magenta'
		if "2014" in label or "2015" in label: m='v'
		if "2016" in label or "2017" in label: m='o'
		plots[label] = ax.scatter(df['x'][i], df['y'][i], c=c, edgecolors='k', marker=m)
	for spine in ax.spines.values(): spine.set_alpha(0.2)
	ax.legend((plots["2014-09"],plots["2015-06"],plots["2016-02"],plots["2017-02"]),
		("2014-09","2015-06","2016-02","2017-02"),
		numpoints=1, loc='upper left', ncol=2, framealpha=1, frameon=True, facecolor='w', columnspacing=1, handlelength=0)
	ax.set_xlabel("PC1 ("+str(int(100*PC1_expl))+"%) variation explained")
	ax.set_ylabel("PC2 ("+str(int(100*PC2_expl))+"%) variation explained")
	ax.set_title("PCoA of Weighted Unifrac matrix (Site A)")


def draw_archaea_percent(data1, data2, ax):
	df1={"date":[], "vals":[]}; df2={"date":[], "vals":[]};
	for year in range(2014,2018):
		for month in range(1,13):
			if year==2014 and month<6: continue
			if year==2017 and month>04: continue
			if month<10: date=str(year)+"-0"+str(month)
			else: date=str(year)+"-"+str(month)
			
			if date in data1:
				for point in data1[date]:
					df1["date"].append(date)
					df1["vals"].append(point)
			if date in data2:
				for point in data2[date]:
					df2["date"].append(date)
					df2["vals"].append(point)
			else:
				df1["date"].append(date); df2["date"].append(date)
				df1["vals"].append(0); df2["vals"].append(0)

	sns.boxplot(x="date", y='vals', data=df1, width=2, linewidth=1, ax=ax, color="magenta")
	sns.boxplot(x="date", y='vals', data=df2, width=2, linewidth=1, ax=ax, color="red")
	
	sns.swarmplot(x="date", y='vals', data=df1, size=3, edgecolor="black", linewidth=1, ax=ax, color="magenta")
	sns.swarmplot(x="date", y='vals', data=df2, size=3, edgecolor="black", linewidth=1, ax=ax, color="red")
	
	ax.set_xticks([3,12,20,25,28,32])
	for tick in ax.get_xticklabels(): tick.set_rotation(45)
	ax.set_ylim(30,100)
	ax.set_ylabel("Relative Archaea abundance")
	ax.set_title("Relative abundance of Archaea in at two sampling sites")
  	ax.axvline(x=15, ls='--')

	legend_elements = [Patch(facecolor='magenta', edgecolor='k', label='Site A'), Patch(facecolor='red', edgecolor='k', label='Site B')]
	ax.legend(handles=legend_elements, loc='lower left', framealpha=1, frameon=True, facecolor='w')
	for spine in ax.spines.values(): spine.set_alpha(0.2)


def draw_clustermap(filename, ax, v, h):
	im = Image.open(filename)
	plt.figimage(im, h, v)
	ax.axis("off")
	ax.set_title("Dissimilarity clustering (Site A)")
	legend_elements = [Patch(facecolor='red', edgecolor='k', label='2014-09'), 
		Patch(facecolor='yellow', edgecolor='k', label='2015-06'),
		Patch(facecolor='cyan', edgecolor='k', label='2016-02'),
		Patch(facecolor='magenta', edgecolor='k', label='2017-02')]
	ax.legend(handles=legend_elements, loc=[-0.05, -0.13], framealpha=1, frameon=True, facecolor='w', ncol=4, columnspacing=0.5, handlelength=1)



##################   START SCRIPT     ######################


# main figure layout:
plt.rc('font', family='sans-serif')
sns.set_palette("husl")
sns.set_style("dark")
fig = plt.figure(figsize=(16, 11))
outer = gridspec.GridSpec(3, 1, wspace=0.15, hspace=0.6, height_ratios=[2,1,1], width_ratios=[1])


for i, label in enumerate(['TOP PANEL','MIDDLE','BOTTOM']):
	if label=="TOP PANEL":
		# make other figures in top row
		inner = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[i], wspace=0.2, width_ratios=[1,0.7,1])
		for j, panel in enumerate(['A', 'B','C']):
			ax = plt.Subplot(fig, inner[j])
			ax.annotate(panel, xy=(-0.1, 1.1), xycoords="axes fraction", fontsize=16)
			if panel=="A":
				print "making PCoA plot..."
				PC1_expl,PC2_expl,PC1,PC2,labels = load_pca_coord("weighted_unifrac_pcoa.txt")
				df = convert_lists_to_df(PC1, PC2, labels)
				draw_pcoa(PC1_expl, PC2_expl, df, ax)
			if panel=="B":
				print "making clustermap plot..."
				os.system("python clustermap.py weighted_unifrac_matrix.tab 4.2")
				vertical=660
 				horisontal=650
				draw_clustermap("clustermap.png", ax, vertical, horisontal)
			if panel=="C":
				print "making Archaea abundance plot..."
				inputfile1="siteA_otu_table.tab"
				inputfile2="siteB_otu_table.tab"
				df1 = load_otu_data(inputfile1, 0)
				df2 = load_otu_data(inputfile2, 0)
				data1=load_taxa_data(df1, "Archaea")
				data2=load_taxa_data(df2, "Archaea")
				draw_archaea_percent(data1, data2, ax)

			fig.add_subplot(ax)

	if label=="MIDDLE" or label=="BOTTOM":
		# make phyla abundance figures
		if label=='MIDDLE': inputfile="siteA_otu_table.tab"
		if label=='BOTTOM': inputfile="siteB_otu_table.tab"
		print "making phyla abundance plot ("+inputfile+")..."
	
		rank_of_interest=2
		taxa_of_interest=["Cyanobacteria","Chloroplast","Cytophagia","Halobacteria"]
	
		df = load_otu_data(inputfile, rank_of_interest)
		taxa_data={}
		for taxa in taxa_of_interest: 
			taxa_data[taxa]=load_taxa_data(df, taxa)

		inner = gridspec.GridSpecFromSubplotSpec(1, len(taxa_of_interest), subplot_spec=outer[i], wspace=0.3)
		alphabet=["A","B","C","D","E","F","G","H","I","J","K","L"]; pos=3
		for j, taxa in enumerate(taxa_data):
			ax = plt.Subplot(fig, inner[j])
			ax.annotate(alphabet[pos], xy=(-0.1, 1.1), xycoords="axes fraction", fontsize=16)
			pos+=1
			data=taxa_data[taxa]
			draw_boxplots(data, taxa, ax, i)
			if j==0 and label=="MIDDLE": ax.set_ylabel("Abundance (%) in Site A")
			if j==0 and label=="BOTTOM": ax.set_ylabel("Abundance (%) in Site B")
			fig.add_subplot(ax)

plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)	
plt.savefig("figure_1.png")
plt.show()



