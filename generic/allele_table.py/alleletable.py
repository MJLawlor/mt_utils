#!/bin/python/

import sys
import os
alleletable=[]
headers = 'CHROM POS REF ALT'.split()

uniquevariants=open(sys.argv[1],'r+')
for line in uniquevariants:
	varline=line.split()
	alleletable.append(varline)
	
def updateTable(file):
	array=[]
 	with open(file) as fh:
		sample = os.path.splitext(file)[0]
		headers.append(sample)
		for line in fh:
			subline=line.split()
			array.append(subline)
	for j in range(len(alleletable)):
		for i in range(len(alleletable)):
			if array[i][1]==alleletable[j][1]:
				if array[i][3]==alleletable[j][3]:
					alleletable[j].append(array[i][3])

def addRef():
	lengths=[]
	for k in range(len(alleletable)):	
		lengths.append(len(alleletable[k]))
        number=max(lengths)
	for l in range(len(alleletable)):	
		if len(alleletable[l])<number:
			alleletable[l].append('REF')

for index in range(2,len(sys.argv)):
	input=sys.argv[index]
	updateTable(input)
	addRef()
		
print ' '.join(headers)

for line in alleletable:
	lengths2=[]
	lengths2.append(len(line))

if max(lengths2)!=min(lengths2):
	raise ValueError("Problem while creating array. Rows not of equal length")

for line in alleletable:
	print (' '.join(line))
