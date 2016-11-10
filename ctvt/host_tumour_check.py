#!/usr/bin/env python
import sys
import os
import re
# Initialise empty arrays
array1, array2, array3, array4 =([] for i in range(4))
def makeArrays(table,array):
        for line in table:
                splitline=line.split()
                array.append(splitline)

file1=sys.argv[1] #file1 should be tumour
tumour=open(file1, 'r+')
file2=sys.argv[2] #file2 should be host
host=open(file2, 'r+')
for line in host:
                splitline=line.split()
                array4.append(splitline[1])
                array2.append(splitline)

makeArrays(tumour, array1)
#makeArrays(host, array2)
for x in range(len(array1)):
        for y in range(len(array2)):
                if array1[x][1]==array2[y][1] and array1[x][2]==array2[y][2] and array1[x][5]>='0.9':
                        array3.append(array1[x])
                elif array1[x][1]==array2[y][1] and array1[x][2]!=array2[y][2]:
                        array3.append(array1[x])
        if array1[x][1] in array4:
                pass
        else:
                array3.append(array1[x])

for line in array3:
        print (' '.join(line))
