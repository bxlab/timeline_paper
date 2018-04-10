#!/usr/bin/env python
import sys

total=[]

for f in sys.argv[1:]:
	a=[]
	for line in open(f): 
		a.append(float(line.strip()))

	if total==[]:
		for i in a: total.append(i)	
	else:
		for i in range(len(a)):
			total[i]+=a[i]

for i in range(len(total)):
	print total[i]/5









