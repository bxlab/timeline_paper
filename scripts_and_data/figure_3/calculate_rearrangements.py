#!/usr/bin/env python

print "importing libraries..."
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


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
	

def write_to_file(filename, list1, list2):
	file = open(filename, 'w') 
	file.write("#abundance\trearangement\n")
 	for i, a in enumerate(list1):
		file.write(str(a)+"\t"+str(list2[i])+"\n")
	file.close()


def find_pathway_carriers(file_name, level):
	print "finding which contigs carry each pathway function..."
	function_carriers={}
	for line in open(file_name):
		if line[0]=="#": continue
		cut=line.strip().split("\t")
		contig = cut[0].split("-")[0]
		if contig not in contigs: continue
		if cut[3]=="NA": continue

		if level==0:
			paths = [cut[1]]
		else:
			paths = cut[4].split("|")
			for i, path in enumerate(paths):
				paths[i] = path.split(";")[-level]

		for path in paths:
			if path not in function_carriers:
				function_carriers[path]=[contig]
			else:
				function_carriers[path].append(contig)
	return function_carriers


def compute_differences_between_two_samples(sample_A, sample_B):
	print "computing representation differences between "+sample_A+" and "+sample_B+"..."
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

	#print "Average rearangement index = " + str(weighted_average(rearangements, total_abund)) + "+/-" + str(weighted_stdev(rearangements, total_abund))
	return functions, total_changes, rearangements, total_abund


# MAIN
depths, samples, contigs = load_contig_depths("contig_depth.tab")
function_carriers = find_pathway_carriers("img_annotation.master", int(sys.argv[1]))

years=["2014", "2015", "2016", "2017"]
for sampleA in years:
	sampleB=str(int(sampleA)+1)
	if sampleB=="2018": sampleB="2014"
	functions, total_changes, rearangements, abundances = compute_differences_between_two_samples(sampleA, sampleB)
	write_to_file(sampleA+"-"+sampleB+".tab", abundances, rearangements)




