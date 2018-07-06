#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np

barplot_data={}
lvl=4

for file in sys.argv[1:]:
	print "loading in "+file
	tally={}
	tot=0
	for line in open(file):
		if line.startswith("NODE"):
			#these are contigs
			continue
		else:
			#these are reads...
			cut=line.strip().split("\t")
			taxa=cut[1].split(";")
			if len(taxa)<lvl+1: continue
			if taxa[1]!="cellular organisms": continue
			tot+=1
			tax = taxa[lvl]
			if tax=="unclassified Archaea (miscellaneous)": tax="Nanohaloarchaea"
	
			if tax not in tally: tally[tax]=0
			tally[tax]+=1
	tmp_groups=[]	
	for k in sorted(tally, key=tally.get, reverse=True):
		tmp_groups.append(k)
		if k not in barplot_data:
			barplot_data[k]=[]
		barplot_data[k].append(100.0*tally[k]/tot)

	
groups=[]
for k in tmp_groups:
	if len(barplot_data[k])<20: continue
	groups.append(k)

groups.append("Other")
barplot_data["Other"]=[]
for i in range(len(sys.argv[1:])):
	tot=0
	for k in groups:
		if k=="Other": continue
		tot+=barplot_data[k][i]
	barplot_data["Other"].append(100.0-tot)


ind = np.arange(20)
width = 0.9
labels = sys.argv[1:]
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
#plt.yscale('log')
plt.yticks(np.arange(0, 101, 10))

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title='Phyla', bbox_to_anchor=(1, 1))

plt.tight_layout(rect=[0,0,0.5,1])

plt.show()

