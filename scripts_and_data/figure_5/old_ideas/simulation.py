#!/usr/bin/env python
import sys
import random
from random import shuffle
from operator import itemgetter


def generate_cell(probabilities, list_of_functions):
	functions = []
	for i, f in enumerate(list_of_functions):
		for i in range(probabilities[i]):
			functions.append(f)
	n = random.randint(2,3)
	cell=[]
	while len(cell)<n:
		func = random.choice(functions)
		if func in cell: continue
		cell.append(func)
	return cell


def generate_community(functions):
	community={}
	bact_abund = [20,14,10,6,3,2,1]
	for abund in bact_abund:
		while True:
			cell = generate_cell([8,6,4,2,1], functions)
			name = "BACT_"+"_".join(cell)
			if name not in community: break
		community[name] = abund

	arch_abund = [10,7,5,3,2,1,1]
	for abund in arch_abund:
		while True:
			cell = generate_cell([1,2,4,6,8], functions)
			name = "ARCH_"+"_".join(cell)
			if name not in community: break
		community[name] = abund
	return community


def calculate_functional_potential(community):
	potential={}
	for cell in community:
		for f in cell.split("_")[1:]:
			if f not in potential:
				potential[f]=0
			potential[f]+=community[cell]
	return potential


def slow_shift(community):
	new_community={}
	for cell in community:
		if cell.startswith("ARCH"):
			new_community[cell] = community[cell]*2
		if cell.startswith("BACT"):
			new_community[cell] = community[cell]/2
	return new_community


def rearrange(community):
	arch_abund=[]
	bact_abund=[]
	for cell in community:
		if cell.startswith("ARCH"):
			arch_abund.append(community[cell])
		if cell.startswith("BACT"):
			bact_abund.append(community[cell])
	shuffle(arch_abund)
	shuffle(bact_abund)
	a=0; b=0
	new_community={}
	for cell in community:
		if cell.startswith("ARCH"):
			new_community[cell] = arch_abund[a]
			a+=1
		if cell.startswith("BACT"):
			new_community[cell] = bact_abund[b]
			b+=1
	return new_community


def check_if_equal(dict1, dict2):
	distance=0
	for key in dict1:
		distance += abs(dict1[key]-dict2[key])
	if distance<5:
		return True
	else:
		return False


def print_community(community):
	for name in sorted(community):
		print str(community[name])+"\t"+name
	print""


def print_functional_potential(potential):
	for f in sorted(potential):
		print str(potential[f])+"\t"+f
	print ""


def print_all_states(community1, community2, community3):
	fields=[]
	for key in sorted(community1):
		if not key.startswith("ARCH"): 
			continue
		fields.append([community1[key], community2[key], community3[key], key])
	for l in sorted(fields, key=itemgetter(0), reverse=True):
		print '\t'.join(map(str,l))
	fields=[]
	for key in sorted(community1):
		if not key.startswith("BACT"):
			continue
		fields.append([community1[key], community2[key], community3[key], key])
	for l in sorted(fields, key=itemgetter(0), reverse=True):
		print '\t'.join(map(str,l))

	
		

functions=['cross', 'diamond', 'droplet', 'star', 'triangle']

i=0
while True:
	start_community = generate_community(functions)
	start_functional_potential = calculate_functional_potential(start_community)

	slow_community = slow_shift(start_community)
	slow_functional_potential = calculate_functional_potential(slow_community)

	fast_community = rearrange(slow_community)
	fast_functional_potential = calculate_functional_potential(fast_community)

	if check_if_equal(fast_functional_potential, slow_functional_potential)==True:
		print str(i)+"\n"

		print_functional_potential(start_functional_potential)
		print_functional_potential(slow_functional_potential)
		print_functional_potential(fast_functional_potential)

		#print_community(start_community)
		#print_community(slow_community)
		#print_community(fast_community)
		print_all_states(start_community, slow_community, fast_community)
		break
	i+=1




