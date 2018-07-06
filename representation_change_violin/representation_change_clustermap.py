#!/usr/bin/env python
import sys
from scipy.stats import pearsonr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

years=["2014", "2015", "2016", "2017"]

reshuff = {}
for sampleA in years:
	reshuff[sampleA]={}
	for sampleB in years:
		if sampleA == sampleB:
			reshuff[sampleA][sampleB]=0
		else:
			f=sampleA+'_'+sampleB+'.tab'
			tot=0
			ct=0
			for line in open(f):
				ct+=1
				tot+=float(line.strip().split("\t")[1])
			reshuff[sampleA][sampleB]=tot/ct


df = pd.DataFrame.from_dict(reshuff)
sns.clustermap(df, col_cluster=True)


plt.show()

