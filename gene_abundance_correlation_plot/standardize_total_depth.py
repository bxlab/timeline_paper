#!/usr/bin/env python
import sys
data=[]
tot=0
ct=0
for line in open(sys.argv[1]):
	n=float(line.strip())
	tot+=n
	ct+=1
	data.append(n)

for i in data:
	print 1000000*i/tot
