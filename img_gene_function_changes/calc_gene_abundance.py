#!/usr/bin/env python
# ./calc_potassium_pump_abundance.py Bacteroidetes "glycine betaine"

import sys
import matplotlib.pyplot as plt
from scipy import stats

# load contig depths
print "loading contig depths"
depths={}
for line in open("contig_depth.tab"):
        cut =line.strip().split("\t")
        if cut[0]=="#contig":
                head = cut
        else:
                depths[cut[0]]={}
                for i in range(1, len(head)):
                        sample=head[i]
                        val=float(cut[i])
                        depths[cut[0]][sample]=val

#import sample read ct
reads={}
for line in open("sample_read_count.tab"):
	if line[0]=="#": continue
	cut=line.strip().split("\t")
	reads[cut[0]]=int(cut[1])

# load taxonomy info
taxa={}
for line in open("contig_taxonomy.tab"):
        cut=line.strip().split("\t")
        #taxa[cut[0]]=cut[1].split(";")
	taxa[cut[0]]=cut[1]

# quant occurances
print "quantifying gene products"
count={}
totals={}
for sample in head[1:]:
	count[sample]=0
	totals[sample]=0



taxa_of_interest=sys.argv[1]
protein_of_interest=sys.argv[2]

for line in open("img_annotation.products"):
	cut =line.strip().split("\t")
	contig=cut[0].split("-")[0]
	if contig not in depths: continue
	if contig not in taxa: continue
	tax=taxa[contig]
	if taxa_of_interest not in tax: continue
	funct=cut[1]
	for sample in head[1:]:
		abund=depths[contig][sample]
		totals[sample]+=abund
		if protein_of_interest in funct:
			count[sample]+=abund

data=[[],[],[],[]]
for sample in sorted(count):
        #abund = count[sample]
	abund = 1000000.0* count[sample] / totals[sample]
	#abund = 1000000.0* count[sample] / reads[sample]
	year_bin=int(sample.split("-")[1])-2014
	data[year_bin].append(abund)

test=stats.ttest_ind(data[1], data[2])
pval=test[1]

plt.style.use('ggplot')
ax = plt.subplot(111)
ax.boxplot(data, patch_artist=True,
	boxprops=dict(facecolor="c", color="k"),
	capprops=dict(color="k"),
	whiskerprops=dict(color="k"),
	flierprops=dict(color="k", markeredgecolor="k"),
	medianprops=dict(color="k"))
#ax.grid(b=True, which='both', color='0.65', linestyle='-')
#ax.set_facecolor("white")
ax.set_xticklabels(['2014', '2015', '2016', '2017'])
ax.set_title(protein_of_interest + ' potential in ' + taxa_of_interest + " (pval= " + str(pval)[:4] + ")")
ax.set_ylabel("Relative abundance of " + protein_of_interest)

plt.show()





