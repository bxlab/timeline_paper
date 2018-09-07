#!/usr/bin/env python
import sys

ct=0


def remove_references(line):
	if "[" in line or "]" in line:
		out=""
		comm=False
		for c in line:
			if c=="[": comm=True
			if comm==False: out+=c
			if c=="]": comm=False
		line=out
	return line


for line in open(sys.argv[1]):
	line=line.strip()
	if line=="": continue
	line = remove_references(line)
	cut = line.split()
	ct += len(cut)


print ct
