#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import operator
from scipy import stats
if len(sys.argv)<2: sys.argv.append("cellular")


def average_iep_statistics(data, labels):
	for i,s1 in enumerate(labels):
		for j,s2 in enumerate(labels):
			if s1==s2: continue
			if int(s1)>int(s2): continue
			test=stats.ttest_ind(data[i], data[j])	
			print s1, s2, test
	
def iep_distribution_statistics(data):
	print "Kolmogorov-Smirnov statistics:"
	for s1 in sorted(data):
		for s2 in sorted(data):
			if s1==s2: continue
			if int(s1.split("-")[1])>int(s2.split("-")[1]): continue
			if int(s1.split("_")[0].split("-")[3])>int(s2.split("_")[0].split("-")[3]): continue
			d1=[]; d2=[]
			for tup in data[s1]:
				iep=tup[0]
				abund=int(tup[1])
				for a in range(abund):
					d1.append(iep)
			for tup in data[s2]:
				iep=tup[0]
				abund=int(tup[1])
				for a in range(abund):
					d2.append(iep)
			
			test=stats.ks_2samp(d1, d2)
			print s1, s2, test


def load_contig_taxonomy(dir):
	print "loading taxonomy..."
	taxonomy={}
	for f in os.listdir(dir):
		sample=f.split('.')[0]
		taxonomy[sample]={}
		for line in open(dir+"/"+f):
			cut=line.strip().split("\t")
			taxonomy[sample]["_".join(cut[0].split("_")[:4])]=cut[1]
	return taxonomy


def load_gene_taxonomy(contig_taxonomy, dir):
	print "loading gene taxonomy..."
	taxonomy={}
	for f in os.listdir(dir):
		sample=f.split('.')[0]
		taxonomy[sample]={}
		for line in open(dir+"/"+f):
			cut=line.strip().split("\t")
			contig = cut[0]
			if contig not in contig_taxonomy[sample]: continue
			taxa = contig_taxonomy[sample][contig]
			gene = cut[8].split(";")[0].split("=")[1]
			taxonomy[sample][gene] = taxa
	return taxonomy


def load_contig_depths(dir):
	print "loading gene depths..."
	depths={}
	for f in os.listdir(dir):
		sample=f.split('.')[0]
		depths[sample]={}
		for line in open(dir+"/"+f):
			cut=line.strip().split("\t")
			contig = "_".join(cut[2].split("_")[:4])
			if int(contig.split("_")[3])<500: continue
			depths[sample][cut[0]] = float(cut[1])
	return depths


def load_gene_ieps(dir, depths, taxonomy):
	print "loading IEP data..."
	ieps={}
	for f in os.listdir(dir):
		sample=f.split('.')[0]
		year=sample.split("-")[1]
		ieps[sample]=[]
		
		for line in open(dir+"/"+f):
			cut=line.strip().split("\t")
			if len(cut)<5 or cut[0]=="Pro_id": continue
			iep=float(cut[-1])
			gene=cut[0]
			if gene not in depths[sample]: continue
			depth = depths[sample][gene]
			if gene in taxonomy[sample]:
				taxa = taxonomy[sample][gene]
			else:
				taxa="cellular"
			ieps[sample].append( (iep, depth, taxa) )
	return ieps


def load_co_taxonomy(tax_file):
	print "\nloading co-assembly contig taxonomy..."
	taxonomy={}
	for line in open(tax_file):
		cut=line.strip().split("\t")
		taxonomy[cut[0]]=cut[1]
	return taxonomy


def load_co_depths(file_name):
	print "loading co-assembly contig depths..."
	depths=np.loadtxt(open(file_name), delimiter="\t", skiprows=0, usecols=range(1,21))
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


def select_contigs_by_taxonomy(contigs, taxonomy, taxa):
	contigs_out={}
	for contig in taxonomy:
		if taxa not in taxonomy[contig]: continue
		if contig not in contigs: continue
		contigs_out[contig]=contigs[contig]
	return contigs_out
	

def get_trk_abundance(annotation_file, depths, samples, contigs):
	print "pulling out Trk abundances..."
	trk_out={}
	depth_totals={}
	for sample in samples:
		trk_out[sample]=0
		depth_totals[sample]=0
	for line in open(annotation_file):
		if line[0]=="#": continue
		cut=line.strip().split("\t")
		contig=cut[0].split("-")[0]
		if contig not in contigs: continue
		for sample in samples:
			depth=depths[contigs[contig]][samples[sample]]
			depth_totals[sample] += depth
			if "potassium uptake protein" in cut[2]:
				trk_out[sample] += depth
	for sample in samples:
		trk_out[sample] = 1000 * trk_out[sample] / depth_totals[sample]
	return trk_out


def get_iep_distribution(ieps):
	print "placing IEPs into distribution lists..."
	isoelectric_point_bins = {}
	for sample in ieps:
		isoelectric_point_bins[sample]=[]
		for i in np.arange(0.0, 14, 0.2):
			isoelectric_point_bins[sample].append(0)
		for tup in ieps[sample]:
			iep = tup[0]
			depth = tup[1]
			iep_bin = int(5*iep)
			isoelectric_point_bins[sample][iep_bin]+=depth
	return isoelectric_point_bins


def group_distributions(iep_dists):
	print "combining IEP distribution lists..."
	summ={}
	for sample in iep_dists:
		new=sample.split("-")[1]
		if new not in summ:
			summ[new]=iep_dists[sample]
		else:
			summ[new]=map(operator.add, summ[new], iep_dists[sample])
	return summ


def merge_distributions(iep_dists):
	print "combining IEP distribution lists..."
	merge=[]
	for sample in iep_dists:
		if merge==[]:
			merge = iep_dists[sample]
		else:
			merge = map(operator.add, merge, iep_dists[sample])
	return merge


def standardize_distributions(iep_dists):
	print "standardizing IEP distributions..."
	std_out={}
	for sample in iep_dists:
		std_out[sample]=[]
		tot=np.sum(iep_dists[sample])
		for i in range(len(iep_dists[sample])):
			std_out[sample].append(100*iep_dists[sample][i]/tot)
	return std_out


def subset_by_taxa(ieps, taxa):
	print "filtering for "+taxa+" genes..."
	subset={}
	for sample in ieps:
		subset[sample]=[]
		for tup in ieps[sample]:
			if taxa in tup[2]:
				subset[sample].append(tup)
	return subset


def weighted_average_iep(list_of_tuples):
	total=0
	weights=0
	for tup in list_of_tuples:
		total += tup[0]*tup[1]
		weights += tup[1]
	return total/weights


def draw_dist(data, ax, colors):
	print "plotting pI curve..."
	for sample in sorted(data):
		if "2014" in sample: c=colors[0]
		elif "2015" in sample: c=colors[1]
		elif "2016" in sample: c=colors[2]
		elif "2017" in sample: c=colors[3]
		elif "Bact" in sample: c="red"
		elif "Halo" in sample: c="blue"
		else: c="green"
		if sample.split("_")[0].split("-")[-1]=='1' or sample.startswith("201") or "SG" not in sample:
			ax.plot(np.arange(0.0, 14, 0.2), data[sample], linewidth=2, color=c, label=sample)
		else:
			ax.plot(np.arange(0.0, 14, 0.2), data[sample], linewidth=2, color=c)
	ax.legend()
	for spine in ax.spines.values(): spine.set_alpha(0.2)
	ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
	ax.set_xlim(3,13)


def draw_iep_boxplot(ieps, ax, c):
	df={"x":[], "y":[]}
	for sample in ieps:
		iep = weighted_average_iep(ieps[sample])
		year = sample.split("-")[1]
		df["x"].append(year)
		df["y"].append(iep)
	sns.boxplot(x="x", y='y', data=df, width=0.5, linewidth=1, ax=ax, palette=c)
	sns.swarmplot(x="x", y='y', data=df, size=5, edgecolor="black", linewidth=0.5, ax=ax, palette=c)
	ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
	for spine in ax.spines.values(): spine.set_alpha(0.2)


def draw_trk_boxplot(trk_dict, ax, c):
	df={"x":[], "y":[]}
	for sample in trk_dict:
		val = trk_dict[sample]
		year = sample.split("-")[1]
		df["x"].append(year)
		df["y"].append(val)
	sns.boxplot(x="x", y='y', data=df, width=0.5, linewidth=1, ax=ax, palette=c)
	sns.swarmplot(x="x", y='y', data=df, size=5, edgecolor="black", linewidth=0.5, ax=ax, palette=c)
	ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
	for spine in ax.spines.values(): spine.set_alpha(0.2)


# START MAIN SCRIPT

# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
sns.set_palette("colorblind")
fig = plt.figure(figsize=(14, 10))
colors=["gold", "cyan", "royalblue", "magenta"]
title_font=16


# loading individual assembly data
contig_taxonomy = load_contig_taxonomy("contig_taxonomy")
gene_depths = load_contig_depths("gene_depths")
gene_taxonomy = load_gene_taxonomy(contig_taxonomy, "contig_annotation")
ieps = load_gene_ieps("gene_ieps", gene_depths, gene_taxonomy)

# loading co-assembly data
co_taxonomy = load_co_taxonomy("co-assembly/contig_taxonomy.tab")
co_depths, co_samples, co_contigs = load_co_depths("co-assembly/contig_depth.tab")


print "\nplotting metagenome IEP curves..."
ax = fig.add_subplot(231)
iep_dists = get_iep_distribution(ieps)
joint_dist = group_distributions(iep_dists)
std_dist = standardize_distributions(joint_dist)
draw_dist(std_dist, ax, colors)
ax.set_xlabel("Isoelectric Point (pI)")
ax.set_ylabel("Relative abundance")
ax.set_title("Protein pI distribution", fontsize=title_font)
ax.annotate("A", xy=(-0.08, 1.02), xycoords="axes fraction", fontsize=title_font)


print "\nplotting average IEP box plots..."
ax = fig.add_subplot(232)
draw_iep_boxplot(ieps, ax, colors)
ax.set_ylabel("Average isoelectric point")
ax.set_title("Average pI", fontsize=title_font)
ax.annotate("B", xy=(-0.08, 1.02), xycoords="axes fraction", fontsize=title_font)


print "\nplotting total Trk potential..."
ax = fig.add_subplot(233)
trk_data = get_trk_abundance("co-assembly/img_annotation.master", co_depths, co_samples, co_contigs)
draw_trk_boxplot(trk_data, ax, colors)
ax.set_ylabel("Average Trk abundance")
ax.set_title("Potassium uptake potential", fontsize=title_font)
ax.annotate("C", xy=(-0.08, 1.02), xycoords="axes fraction", fontsize=title_font)


print "\nplotting taxa IEP curves..."
ax = fig.add_subplot(234)
dist={}
for i,taxa in enumerate(["Bacteroidetes", "Halobacteria"]):
	taxa_ieps = subset_by_taxa(ieps, taxa)
	iep_dists = get_iep_distribution(taxa_ieps)
	dist[taxa] = merge_distributions(iep_dists)
	std_dist = standardize_distributions(dist)
draw_dist(std_dist, ax, colors)
ax.set_xlabel("Isoelectric Point (pI)")
ax.set_ylabel("Relative abundance")
ax.set_title("Protein pI distribution by taxonomy", fontsize=title_font)
ax.annotate("D", xy=(-0.08, 1.02), xycoords="axes fraction", fontsize=title_font)


print "\nplotting average IEP box plots for Halobacteria..."
ax = fig.add_subplot(235)
taxa_ieps = subset_by_taxa(ieps, "Halobacteria")
draw_iep_boxplot(taxa_ieps, ax, colors)
ax.set_ylabel("Average isoelectric point")
ax.set_title("Average pI of Halobacteria", fontsize=title_font)
ax.annotate("E", xy=(-0.08, 1.02), xycoords="axes fraction", fontsize=title_font)


print "\nplotting total Trk potential of Halobacteria..."
ax = fig.add_subplot(236)
halo_contigs = select_contigs_by_taxonomy(co_contigs, co_taxonomy, "Halobacteria")
trk_data = get_trk_abundance("co-assembly/img_annotation.master", co_depths, co_samples, halo_contigs)
draw_trk_boxplot(trk_data, ax, colors)
ax.set_ylabel("Average Trk abundance")
ax.set_title("Potassium uptake in Halobacteria", fontsize=title_font)
ax.annotate("F", xy=(-0.08, 1.02), xycoords="axes fraction", fontsize=title_font)


plt.tight_layout()
plt.savefig("figure_S5.png", dpi=300)












