#!/usr/bin/env python
import sys
data={}
for line in open(sys.argv[1]):
	cut=line.strip().split("\t")
	if cut[0]=="Category": print line.strip()
	else:
		for i in range(1,len(cut)):
			cut[i]=float(cut[i])
		func=cut[0].split(";")[1]
		if func not in data: data[func]=cut[1:]
		else:
			data[func]=map(sum, zip(data[func],cut[1:]))


for k in data:
	for i in range(len(data[k])):
		data[k][i]=str(data[k][i])
	print k+"\t"+"\t".join(data[k])	
		
