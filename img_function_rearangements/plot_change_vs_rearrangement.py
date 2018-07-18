#!/usr/bin/env python

print "importing libraries..."
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns


def scaled_fold(a, b):
	if a+b == 0:
		return 0
	else:
		return (b-a)/(a+b)


def weighted_average(numbers, weights):
	total=0
	weight_total=0
	for i in range(len(numbers)):
		weight_total+=weights[i]
		total+=numbers[i]*weights[i]
	if weight_total==0 or total==0:
		return 0
	else:
		return total/weight_total


def weighted_stdev(values, weights):
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return math.sqrt(variance)


def artificial_std(values, weights):
	maximum=max(weights)
	out = []
	for i, val in enumerate(values):
		weight = int(1000*weights[i]/maximum)
		for j in range(weight):
			out.append(val)
	return out
		

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
	print "standardizing contig depths in each sample"
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

def write_to_file(filename, list1, list2):
	file = open(filename, 'w') 
	file.write("#abundance\trearangement\n")
 	for i, a in enumerate(list1):
		file.write(str(a)+"\t"+str(list2[i])+"\n")
	file.close()

def find_gene_product_carriers(file_name):
	print "finding which contigs carry each gene product..."
	function_carriers={}
	for line in open(file_name):
		if line[0]=="#": continue
		cut=line.strip().split("\t")
		if cut[0].split("-")[0] not in contigs: continue
		if cut[1] not in function_carriers:
			function_carriers[cut[1]]=[cut[0].split("-")[0]]
		else:
			function_carriers[cut[1]].append(cut[0].split("-")[0])
	return function_carriers

def find_pathway_carriers(file_name):
	print "finding which contigs carry each pathway function..."
	function_carriers={}
	for line in open(file_name):
		if line[0]=="#": continue
		cut=line.strip().split("\t")
		if cut[0].split("-")[0] not in contigs: continue
		if cut[3]=="NA": continue

		paths=cut[3].split(";")
		for path in paths:
			if path not in function_carriers:
				function_carriers[path]=[cut[0].split("-")[0]]
			else:
				function_carriers[path].append(cut[0].split("-")[0])
	return function_carriers


def compute_differences_between_two_samples(sample_A, sample_B):
	if len(sys.argv)>2: print "computing representation differences between "+sample_A+" and "+sample_B+"..."
	functions=[]
	total_changes=[]
	rearangements=[]
	total_abund=[]

	for function in function_carriers:
		if len(function_carriers[function])<2: continue

		total_abundance_A = 0
		total_abundance_B = 0
		indiv_change = []		
		contig_weights=[]

		for contig in function_carriers[function]:
			if len(sample_A.split("-"))>2:
				abund_A = depths[contigs[contig]][samples[sample_A]]
				abund_B = depths[contigs[contig]][samples[sample_B]]
			else:
				abund_A = 0
				abund_B = 0
				for sample in samples:
					if sample_A in sample: abund_A += depths[contigs[contig]][samples[sample]]
					if sample_B in sample: abund_B += depths[contigs[contig]][samples[sample]]

			delta = scaled_fold(abund_A, abund_B)
			indiv_change.append(np.absolute(delta))
			contig_weights.append(abund_A+abund_B)

			total_abundance_A += abund_A
			total_abundance_B += abund_B
	
		total_change = scaled_fold(total_abundance_A, total_abundance_B)
		#rearangement = np.mean(indiv_change)
		rearangement = weighted_average(indiv_change, contig_weights)
		#rearangement = weighted_average(indiv_change, contig_weights) - np.absolute(total_change)

		functions.append(function)
		total_changes.append(total_change)
		rearangements.append(rearangement)
		total_abund.append(total_abundance_A + total_abundance_B)
	
	for i in range(len(total_abund)):
		total_abund[i]=50*total_abund[i]/max(total_abund)

	print "Average rearangement index = " + str(weighted_average(rearangements, total_abund)) + "+/-" + str(weighted_stdev(rearangements, total_abund))
	
	return functions, total_changes, rearangements, total_abund


# MAIN
depths, samples, contigs = load_contig_depths("contig_depth.tab")
depths = standardize_contig_depths(depths, samples, contigs)
function_carriers = find_gene_product_carriers("img_annotation.master")
#function_carriers = find_pathway_carriers("img_annotation.master")

if len(sys.argv)>1:
	sampleA=sys.argv[1]
	sampleB=sys.argv[2]
	functions, total_changes, rearangements, abundances = compute_differences_between_two_samples(sampleA, sampleB)
	
	print "plotting..."
	#setting up subplots
	fig = plt.figure(figsize=(8, 8))
	grid = plt.GridSpec(6, 4, hspace=0, wspace=0)
	plt.title(sampleA+" -> "+sampleB, fontsize = 20)
	plt.axis('off')
	
	# plot scatter
	main_ax = fig.add_subplot(grid[1:, :-1])
	main_ax.scatter(total_changes, rearangements, alpha=0.2, s=abundances, c='k')
	main_ax.set_xlabel("Fold change in total abundance")
	main_ax.set_ylabel("Rearangement index of gene functions")
	main_ax.set_xlim(-1,1)
	main_ax.set_ylim(0,1)
	for xmaj in main_ax.xaxis.get_majorticklocs(): main_ax.axvline(x=xmaj, ls='--', alpha=0.3, c='k', linewidth=0.5)
	for ymaj in main_ax.yaxis.get_majorticklocs(): main_ax.axhline(y=ymaj, ls='--', alpha=0.3, c='k', linewidth=0.5)

	# plot x-axis violin
	x_vio = fig.add_subplot(grid[0, :-1], xticklabels=[], yticklabels=[])
	x_vio = sns.violinplot(x=artificial_std(total_changes, abundances))
	x_vio.set_xlim(-1,1)
	for xmaj in x_vio.xaxis.get_majorticklocs(): x_vio.axvline(x=xmaj, ls='--', alpha=0.3, c='k', linewidth=0.5)

	# plot y-axis violin
	y_vio = fig.add_subplot(grid[1:, -1], xticklabels=[], yticklabels=[])
	write_to_file(sys.argv[1]+"-"+sys.argv[2]+".tab", abundances, rearangements)
	y_vio = sns.violinplot(y=artificial_std(rearangements, abundances))
	y_vio.set_ylim(0,1)
	for ymaj in y_vio.yaxis.get_majorticklocs(): y_vio.axhline(y=ymaj, ls='--', alpha=0.3, c='k', linewidth=0.5)


	plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95)
	#plt.savefig(sampleA+"_vs_"+sampleB+".png")
	plt.show()
	quit()



for sampleA in sorted(samples):
	for sampleB in sorted(samples):
		if sampleA==sampleB:
			print "\t".join([sampleA, sampleB, "0.0", "0.0"])
			continue
		functions, total_changes, rearangements, abundances = compute_differences_between_two_samples(sampleA, sampleB)
		print "\t".join([sampleA, sampleB, str(weighted_stdev(total_changes, abundances)), str(weighted_average(rearangements, abundances))])






