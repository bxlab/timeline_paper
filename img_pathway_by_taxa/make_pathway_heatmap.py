#!/usr/bin/env python

print "importing libraries..."
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns
import pandas as pd


def load_taxonomy(tax_file):
	print "loading contig taxonomy..."
	taxonomy={}
	for line in open(tax_file):
		cut=line.strip().split("\t")
		taxonomy[cut[0]]=cut[1]
	return taxonomy


def load_contig_depths(file_name):
	print "loading contig depths..."
	depths=np.loadtxt(open(file_name), delimiter="\t", skiprows=0, usecols=range(1,21))

	print "loading contig and sample name mappings..."
	samples={}
	contigs={}
	for i, line in enumerate(open("contig_depth.tab")):
		cut=line.strip().split("\t")
		if i==0:
			header=cut
			for i in range(len(header)):
				if header[i][0]=="#": continue
				samples[header[i]]=i-1
		else:
			contigs[cut[0]]=i-1
	
	# using this format: depths[contigs[contig]][samples[sample]
	return depths, samples, contigs
	

def standardize_contig_depths(depths, samples, contigs):
	print "standardizing contig depths in each sample..."
	for sample in samples:
		if len(sys.argv)>2:
			if sys.argv[1] not in sample and sys.argv[2] not in sample:
				continue
		total=0
		for contig in contigs:
			length=int(contig.split("_")[3])
			val = depths[contigs[contig]][samples[sample]] * length
			total += val
		
		depths[:,samples[sample]] /= total
		depths[:,samples[sample]] *=1000000

	return depths


def find_pathway_carriers(file_name):
	print "finding which contigs carry each pathway function..."
	function_carriers={}
	for line in open(file_name):
		if line[0]=="#": continue
		cut=line.strip().split("\t")
		if cut[0].split("-")[0] not in contigs: continue
		if cut[3]=="NA": continue

		paths=cut[4].split("|")
		for path in paths:
			cut_path=path.split(";")
			#if cut_path[0]=="Metabolism": continue

			#path=";".join(cut_path[:2])
			#path=";".join(cut_path[:3])
			path=";".join(cut_path[1:3])

			if path not in function_carriers:
				function_carriers[path]=[cut[0].split("-")[0]]
			else:
				function_carriers[path].append(cut[0].split("-")[0])
	return function_carriers


def get_pathway_abundances(function_carriers, depths):
	print "making total pathway abundance table..."
	pathway_abundances={}
	for sample in samples:
		pathway_abundances[sample]={}
	for function in function_carriers:
		for sample in samples:
			pathway_abundances[sample][function]=0
		for contig in function_carriers[function]:
			for sample in samples:
				depth=depths[contigs[contig]][samples[sample]]
				pathway_abundances[sample][function] += depth
	df = pd.DataFrame.from_dict(pathway_abundances)
	return df


def remove_stoic_pathways(pathway_abundances):
	print "selecting differential abundant pathways..."
	df=pathway_abundances
	# standardize columns by total sum in each column
	df = df.div(df.sum(axis=0), axis=1); df=1000000*df

	tot=0
	for gene in df.index.values:
		t2014=[];t2015=[];t2016=[];t2017=[]
		for sample in df.columns.values:
			if "2014" in sample: t2014.append(df[sample][gene])
			if "2015" in sample: t2015.append(df[sample][gene])
			if "2016" in sample: t2016.append(df[sample][gene])
			if "2017" in sample: t2017.append(df[sample][gene])
		anova = stats.f_oneway(t2014, t2015, t2016, t2017)
		if anova[1]<0.01:
			tot+=1
			print gene
		else:
			df=df.drop(gene)
	print str(tot) + " signifficant functional groups found!"
	return df


def plot_even_clustermap(df, log):
	print "plotting clustermap..."
	lut=[]
	for sample in df.columns.values:
		if "2014" in sample: lut.append('r')
		if "2015" in sample: lut.append('g')
		if "2016" in sample: lut.append('b')
		if "2017" in sample: lut.append('c')

	if log==False:
		df = df.div(df.max(axis=1), axis=0)
		g = sns.clustermap(df, col_cluster=True, col_colors=lut, cmap="Blues", figsize=(12,12), vmin=0, vmax=1)
	if log==True:
		df=df+0.01
		df = np.log(df)
		g = sns.clustermap(df, col_cluster=True, col_colors=lut, cmap="Blues", figsize=(12,12))

	plt.subplots_adjust(left=0.05, bottom=0.1, right=0.6, top=0.95)
	plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
	plt.show()


def select_taxonomy_contigs (function_carriers, taxonomy, taxa):
	print "sorting for "+taxa+" contigs..."
	for function in function_carriers:
		good_contigs=[]
		for contig in function_carriers[function]:
			if contig not in taxonomy: continue
			if taxa in taxonomy[contig]:
				good_contigs.append(contig)
		function_carriers[function]=good_contigs
	return function_carriers

def remove_low_abundance_pathways(df):
	df=df.T
	df.drop([col for col, val in df.sum().iteritems() if val < 1], axis=1, inplace=True)
	return df.T


# MAIN
taxonomy = load_taxonomy("contig_taxonomy.tab")
depths, samples, contigs = load_contig_depths("contig_depth.tab")
depths = standardize_contig_depths(depths, samples, contigs)
function_carriers = find_pathway_carriers("img_annotation.master")
if len(sys.argv)>1: function_carriers = select_taxonomy_contigs(function_carriers, taxonomy, sys.argv[1])

all_pathway_abundances = get_pathway_abundances(function_carriers, depths)
pathway_abundances = remove_low_abundance_pathways(all_pathway_abundances)
pathway_abundances.to_csv("pathway_abundance.tab", sep='\t')

differential_pathways = remove_stoic_pathways(pathway_abundances)

differential_pathways.to_csv("differential_pathways.tab", sep='\t')
plot_even_clustermap(differential_pathways, False)



