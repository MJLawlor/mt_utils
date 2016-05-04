#!/usr/bin/python
# author:       Máire Ní Leathlobhair
# date:         May 2016
# description:  Filters tumour variants against matched normal based on variant allele fraction (VAF). If a variant occurs in both tumour
#               and matched normal with VAF>=0.9, variant is not filtered. If a variant occurs in both tumour and matched normal with
#               VAF<0.9 in tumour, variant is discarded. If a variant is called only in the tumour, variant is not filtered.


import sys
import os
		
file1=open(sys.argv[1], 'r+') ##tumour variant list
file2=open(sys.argv[2], 'r+') ##host variant list

def vafCompare(tumour,host):
   tumourArray=[]
   outputArray=[]
   hostVariants=[]
   hostArray=[]
   for line in host:
        splitline=line.split()
        hostVariants.append(splitline[1])
        hostArray.append(splitline)
   for line in tumour:
        splitline=line.split()
        tumourArray.append(splitline)
    
   for x in range(len(tumourArray)):
        for y in range(len(hostArray)):
		        if tumourArray[x][1]==hostArray[y][1] and tumourArray[x][2]==hostArray[y][2] and tumourArray[x][3]>='0.9':
			          outputArray.append(tumourArray[x])
		        elif tumourArray[x][1]==hostArray[y][1] and tumourArray[x][2]!=hostArray[y][2]:
			          outputArray.append(tumourArray[x])
        if tumourArray[x][1] in hostVariants:
            pass
        else:
            outputArray.append(tumourArray[x])
   for line in outputArray:
       print (' '.join(line))

vafCompare(file1,file2)
