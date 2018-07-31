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
	for i, line in enumerate(open(file_name)):
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
			#print gene
		else:
			df=df.drop(gene)
	print str(tot) + " signifficant functional groups found!"
	return df


def plot_even_clustermap(df):
	sns.set(font_scale=0.6)
	print "plotting clustermap..."
	lut=[]
	for sample in df.columns.values:
		if "2014" in sample: lut.append(sys.argv[3])
		if "2015" in sample: lut.append(sys.argv[4])
		if "2016" in sample: lut.append(sys.argv[5])
		if "2017" in sample: lut.append(sys.argv[6])

	df = df.div(df.max(axis=1), axis=0)
	g = sns.clustermap(df, cmap="Blues", figsize=(float(sys.argv[1]), float(sys.argv[2])), vmin=0, vmax=1,
	yticklabels=False, xticklabels=False, col_cluster=True, col_colors=lut)

	plt.subplots_adjust(left=0.05, bottom=0.1, right=0.6, top=0.95)
	plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)


def remove_low_abundance_pathways(data):
        df=data.copy().T
        df.drop([col for col, val in df.sum().iteritems() if val < 1], axis=1, inplace=True)
        return df.T


def select_taxonomy_contigs (function_carriers, taxonomy, taxa):
        functions_out={}
        print "sorting for "+taxa+" contigs..."
        for function in function_carriers:
                good_contigs=[]
                for contig in function_carriers[function]:
                        if contig not in taxonomy: continue
                        if taxa in taxonomy[contig]:
                                good_contigs.append(contig)
                functions_out[function]=good_contigs
        return functions_out



# MAIN

# prepare data
taxonomy = load_taxonomy("contig_taxonomy.tab")
depths, samples, contigs = load_contig_depths("contig_depth.tab")
function_carriers_ALL = find_pathway_carriers("img_annotation.master")

for taxa in ["cellular", "Bacteria", "Archaea"]:
	print "\nPlotting "+ taxa
	function_carriers = select_taxonomy_contigs(function_carriers_ALL, taxonomy, taxa)
	pathway_abundances = get_pathway_abundances(function_carriers, depths)
	pathway_abundances = remove_low_abundance_pathways(pathway_abundances)
	differential_pathways = remove_stoic_pathways(pathway_abundances)
	plot_even_clustermap(differential_pathways)
	if taxa=="cellular": taxa="All"
	plt.savefig(taxa+".png", bbox_inches="tight", dpi=600)
	plt.clf()










