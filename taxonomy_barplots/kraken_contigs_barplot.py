#!/usr/bin/env python
# ./kraken_contigs_barplot.py assembly.kraken contig_abundance.tab
import sys
import matplotlib.pyplot as plt
import numpy as np

taxonomy={}
for line in open(sys.argv[1]):
	#these are reads...
	cut=line.strip().split("\t")
	taxa=cut[1].split(";")
	taxonomy[cut[0]]=taxa


barplot_data={}
lvl=4

for line in open(sys.argv[2]):
	cut=line.strip().split("\t")
	if line[0]=="#": 
		head=cut
		continue
	if cut[0] not in taxonomy: continue
	taxa=taxonomy[cut[0]]
	if len(taxa)<lvl+1: continue
	if taxa[1]!="cellular organisms": continue
	tax=taxa[lvl]
	if tax=="unclassified Archaea (miscellaneous)": tax="Nanohaloarchaea"
	l=int(cut[0].split("_")[3])
	if tax not in barplot_data: 
		barplot_data[tax]=[]
		for i in range(20): barplot_data[tax].append(0)

	for i in range(1, len(cut)):
		abund=float(cut[i])*l
		barplot_data[tax][i-1]+=abund

eg={}
for k in barplot_data: eg[k] = barplot_data[k][0]
groups=[]
for k in sorted(eg, key=eg.get, reverse=True): groups.append(k)
relevant=groups[:15]
print relevant

barplot_data["Other"]=[]
for i in range(20): barplot_data["Other"].append(0)
for i in range(20):
	tot=0
	for k in barplot_data:
		tot+=barplot_data[k][i]
	for k in barplot_data:
		percent = 100.0*barplot_data[k][i]/tot
		if k in relevant:
			barplot_data[k][i]=percent
		else:
			barplot_data["Other"][i]+=percent
			barplot_data[k][i]=0		

	
relevant.append("Other")
groups=relevant

ind = np.arange(20)
width = 0.9
labels = head[1:]
for i in range(len(labels)):
	label=labels[i]
	labels[i]="-".join(label.split("/")[-1].split(".")[0].split("-")[3])


cum_height=[]
for i in range(len(labels)): cum_height.append(0)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1)


for group in groups:
	b = ax.bar(ind, barplot_data[group], width, bottom=cum_height, label=group)
	for i in range(len(cum_height)):
		cum_height[i]+=barplot_data[group][i]

plt.ylabel('Relative abundance')
plt.title('Community composition across timeline replicates')
plt.xticks(ind, labels)
#plt.yticks(np.arange(0, 101, 10))

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title='Phyla', bbox_to_anchor=(1, 1))

plt.tight_layout(rect=[0,0,0.5,1])

plt.show()

