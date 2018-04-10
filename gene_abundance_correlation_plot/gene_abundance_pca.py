#!/usr/bin/env python
import sys
from scipy.stats import pearsonr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

data={"sample":[]}
for f1 in sys.argv[1:]:
	sample=f1.split(".")[0].split("-")[1]+"_"+f1.split(".")[0].split("-")[3]
	data["sample"].append(sample)
	ct=0
	for line in open(f1):
		gene="gene_"+str(ct)
		if gene not in data: 
			data[gene]=[]
		data[gene].append(float(line.strip()))
		ct+=1
df = pd.DataFrame.from_dict(data)
df = df.set_index('sample')
print "finished loading data into df"

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(df)
principalDf = pd.DataFrame(data = principalComponents, columns = ['component_1', 'component_2'])

df = df.reset_index()
print df
finalDf = pd.concat([principalDf, df[["sample"]]], axis = 1)


var = pca.explained_variance_ratio_

fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1 (variance explained = '+ str(var[0]*100)[:4]+"%", fontsize = 10)
ax.set_ylabel('Principal Component 2 (variance explained = '+ str(var[1]*100)[:4]+"%", fontsize = 10)

colors = {"2014": 'r', "2015":'g', "2016":'b', "2017":'c'}
for i in finalDf.index:
	x = finalDf["component_1"][i]
	y = finalDf["component_2"][i]
	sample = finalDf["sample"][i].split("_")[0]
	color=colors[sample]
	if finalDf["sample"][i].split("_")[1]=="1": ax.scatter(x, y, c=color, s=50, label=sample)	
	else: ax.scatter(x, y, c=color, s=50)

ax.legend()
ax.grid()


plt.show()






