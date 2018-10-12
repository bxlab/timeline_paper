#!/usr/bin/env python
# This script takes in the IMG product KO ID annotation file, a table linking the KO identifiers to the KEGG brite pathways, and a table linking each pathway to its proper name.
# Produces a master file linking each gene ID to its possible pathways. This file (img_annotation.master) is used in subsequent analysis.

import sys

# load the BRITE pathway names
brite2func={}
for line in open(sys.argv[1]):
        cut=line.strip().split("\t")
        if "Drug Dev" in cut[1] or "Human Dis" in cut[1] or "Organismal Sys" in cut[1]: continue
        brite2func[cut[0]]=cut[1]

# load the BRITE pathways that each KO ID can contribute to
ko2brite={}
for line in open(sys.argv[2]):
	cut=line.strip().split("\t")
	brites=cut[1].split(";")
	good_brites=[]
	for i in range(len(brites)):
		brite=brites[i]
		if brite[2:] in brite2func:
			good_brites.append(brite)
	cut[1]=";".join(good_brites)
	if cut[1]!="":
		ko2brite[cut[0]]=cut[1]


# Go over the annotage genes and pring out the KO, BRITE ID, and KEGG pathway name
print "\t".join(["#Locus", "KEGG ID", "Gene Product", "BRITE Pathways", "Functions"])
for line in open(sys.argv[3]):
	cut=line.strip().split("\t")
	print cut
	if "KO" not in cut[2]:
		cut.append("NA")
		cut.append("NA")
	elif cut[2].split(":")[1] in ko2brite:
		cut.append(ko2brite[cut[2].split(":")[1]])
		brites = ko2brite[cut[2].split(":")[1]].split(";")
		functions=[]
		for brite in brites:
			functions.append(brite2func[brite[2:]])
		cut.append("|".join(functions))
	else:
		cut.append("NA")
		cut.append("NA")
	ID=cut[2]; name=cut[1]
	cut[1]=ID; cut[2]=name
	print "\t".join(cut)



