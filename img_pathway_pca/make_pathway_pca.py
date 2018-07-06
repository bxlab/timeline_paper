#!/usr/bin/env python
import sys
from scipy.stats import pearsonr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

df = pd.read_csv(sys.argv[1], delimiter="\t", index_col="Category").T
data = df.div(df.sum(axis=0), axis=1); data=1000000*data


standard = StandardScaler()
Data_SArray = standard.fit_transform(data)
Data_Standard = pd.DataFrame(Data_SArray)

samples=list(data.index.values)
indeces=list(Data_Standard.index.values)
mapper = dict(zip(indeces, samples))
df=Data_Standard.rename(index=mapper)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(df)
principalDf = pd.DataFrame(data = principalComponents, columns = ['component_1', 'component_2'])

df = df.reset_index()
finalDf = pd.concat([principalDf, df[["index"]]], axis = 1)


var = pca.explained_variance_ratio_

fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1 (variance explained = '+ str(var[0]*100)[:4]+"%", fontsize = 10)
ax.set_ylabel('Principal Component 2 (variance explained = '+ str(var[1]*100)[:4]+"%", fontsize = 10)
ax.grid()

colors = {"2014": 'r', "2015":'g', "2016":'b', "2017":'c'}
for i in finalDf.index:
	x = finalDf["component_1"][i]
	y = finalDf["component_2"][i]
	sample = finalDf["index"][i].split("-")[1]
	color=colors[sample]
	if finalDf["index"][i].split("-")[3]=="1": ax.scatter(x, y, c=color, s=50, label=sample)	
	else: ax.scatter(x, y, c=color, s=50)

ax.legend()


plt.show()






