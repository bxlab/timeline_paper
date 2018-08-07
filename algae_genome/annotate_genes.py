#!/usr/bin/env python
# ./annotate_genes.py algae_genome.blast algae_genome.gff
import sys

# import gene names
functions={}
for line in open("uniprot_sprot.fasta"):
	if line[0]==">":
		cut=line[1:-1].split()
		name=cut[0]
		func=""
		for f in cut[1:]:
			if f.startswith("OS="): break
			func += " "+f
		functions[name] = func.strip()


# go through blast file
genes={}
for line in open(sys.argv[1]):
	cut=line.strip().split("\t")
	gene = cut[0].split("|")[0]
	identity = float(cut[2])
	coverage = 100.0*float(cut[3])/float(cut[0].split("|")[2].split("_")[0])
	bit = float(cut[11])

	if identity>40 and coverage>70:
		if gene in genes:
			if genes[gene][3]>bit:
				continue
		genes[gene] = [ functions[cut[1]], identity, coverage, bit]
		


#fix gff file
for line in open(sys.argv[2]):
	if line[0]=="#":
		print line.strip()
		continue
	cut = line.strip().split("\t")
	info = cut[8].split(", ")
	gene_id = "gene_" + str(info[0].split("=")[-1])
	if gene_id in genes:
		product_info = genes[gene_id]
		info.append("product="+str(product_info[0]))
		info.append("aa_identity="+str(product_info[1]))
		info.append("aa_query_cover="+str(product_info[2]))
	else:
		info.append("product=unknown")
		info.append("aa_identity=na")
		info.append("aa_query_cover=na")
		
	cut[8]=", ".join(info)
	print "\t".join(cut)
