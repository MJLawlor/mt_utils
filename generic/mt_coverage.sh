#!/bin/bash

# author:       Máire Ní Leathlobhair (ml677@cam.ac.uk)
# date:         May 2016
# simple script for extracting read coverage from a bam file
# input requited: provide list of bed files in a text file and a bed file with coordinates of region where read coverage should be measured

input="$1" # list of bam files 
bed="$2" # bed file containing MT coordinates
length="$3" # length of MT genome
output="$4" # path to folder where we wish to direct output

bams_coverage=($(awk '{print $1}' $input))
for i in "${bams_coverage[@]}";do
	echo $i
	echo -n "$i " >> ${output}coverage.dat
	samtools depth -q 0 -Q 0 -b $bed ${i}.bam > ${output}${i}.dat
 	awk -v l=$length '{sum+=$3; $4=sum/l}{print $4}' ${output}${i}.dat | tail -n 1 >> ${output}coverage.dat
done
