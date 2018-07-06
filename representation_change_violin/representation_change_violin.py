#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd


sns.set_style("whitegrid")
samples=["2014", "2015", "2016", "2017"]

data={}
for sampleA in samples:
	sampleB=str(int(sampleA)+1)
	if sampleA=="2017": sampleB="2014"
	#if int(sampleA)>int(sampleB): continue
	f=sampleA+"_"+sampleB+".tab"
	A=[]
	B=[]
	for line in open(f):
		cut=line.strip().split("\t")
		a=float(cut[0])
		b=float(cut[1])
		A.append(a)
		B.append(b)
	#axarr[x,y].set_yscale('log')
	data[sampleA+"->"+sampleB] = B



		#axarr[x,y].plot(A,B, 'ko', alpha=0.05, ms=1)
		#if x==0: axarr[x,y].set_title(sampleB)
		#if y!=0: axarr[x,y].get_yaxis().set_visible(False)
		#if x!=3: axarr[x,y].get_xaxis().set_visible(False)
		#if y==0: axarr[x,y].set_ylabel(sampleA, fontsize=14)

		#axarr[x,y].set_ylim(0, 1)
		#axarr[x,y].set_xlim(-1, 1)
data_df = pd.DataFrame.from_dict(data)
ax = sns.violinplot(data=data_df)

#plt.savefig("figure.png")
#plt.tight_layout(pad=0.1, w_pad=0.1, h_pad=0.1)
plt.show()
