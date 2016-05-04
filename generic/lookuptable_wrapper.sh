#!/bin/bash

# mitotypus data post-processing pipeline
# author:       Máire Ní Leathlobhair (ml677@cam.ac.uk)
# date:         May 2016

# input
# $1: absolute path to folder containing Mitotypus output
# $2: total number of samples (including host and tumour)
# $3: absolute path to scripts folder
# $4: absolute path to MT reference genome

if [ "$#" -ne 4 ]; then
	echo "| lookuptable_wrapper.sh: Creates variant lists and VAF plots from Mitotypus output"
	echo "|"
	echo "| Usage: lookuptable_wrapper.sh /path/to/GenotypedVariants_final.vcf [number of samples to be analysed] /path/to/scripts/folder /path/to/reference/MT/genome"
fi

PLATYPUSOUT="$1"
PIPELINE="${1}/pipeline"
VARIANTS="${1}/pipeline/variant_lists"
SAMPLES="$2"
SCRIPTS="$3"
LOGS="${1}/pipeline/logs"
DATE=$(date +%y-%m-%d)
REFERENCE="$4"
mkdir -p ${PIPELINE} ${VARIANTS} ${VARIANTS}/tumour ${VARIANTS}/host ${VARIANTS}/tumour/substitutions ${VARIANTS}/host/substitutions ${VARIANTS}/host/indels ${VARIANTS}/tumour/indels ${VARIANTS}/tumour/substitutions/pre-lookup ${VARIANTS}/tumour/substitutions/post-lookup ${PIPELINE}/vaf 

function makeUnfilteredVariantLists { 
	echo -e "\n(1) EXTRACTING VARIANT LISTS FROM GenotypedVariants_final.vcf\n"
	input=${PLATYPUSOUT}/GenotypedVariants_final.vcf
	loopvar=$(($SAMPLES + 9))
#  make temp vcfs containing header and variant call data
  egrep -v '#' "$input" > ${PLATYPUSOUT}/data.tmp.vcf
	egrep -v '##' "$input" | head -1  > ${PLATYPUSOUT}/sample_names.tmp.vcf
	for i in $(seq 10 "$loopvar");
	do 
		name=$(awk -v y="$i" '{print $y}' ${PLATYPUSOUT}/sample_names.tmp.vcf)
		echo -e "\nGenerating variant list for $name\n"	
		awk -v z="$i" '{print $1,$2,$4,$5,$6,$7,$z}' ${PLATYPUSOUT}/data.tmp.vcf |sed -e 's/:/ /g'| awk '{if($12!="0"){$13=$12 / $11 ;print $1,$2,$3">"$4,$5,$6,$13}}' | awk '{if(length($3) == 3){print 
$0}}' > ${VARIANTS}/${name}_substitutions.txt
		awk -v z="$i" '{print $1,$2,$4,$5,$6,$7,$z}' ${PLATYPUSOUT}/data.tmp.vcf |sed -e 's/:/ /g'| awk '{if($12!="0"){$13=$12 / $11 ;print $1,$2,$3">"$4,$5,$6,$13}}' | awk '{if(length($3) > 3){print $
0}}' > ${VARIANTS}/${name}_indels.txt
	done
	## sort variant lists into folders
	mv ${VARIANTS}/*H*_indels.txt ${VARIANTS}/host/indels
	mv ${VARIANTS}/*H*_substitutions.txt ${VARIANTS}/host/substitutions/
	mv ${VARIANTS}/*T*_indels.txt ${VARIANTS}/tumour/indels
	mv ${VARIANTS}/*T*_substitutions.txt ${VARIANTS}/tumour/substitutions/pre-lookup
	###clean up
	rm ${PLATYPUSOUT}/*tmp*
}

function lookupStep {
# compare variants called in matched host tumour samples by PM value and discard host contamination
 echo -e "\n(2) FILTERING VARIANT LISTS. EXECUTING LOOKUP STEP\n"
	
	for x in $(ls ${VARIANTS}/tumour/substitutions/pre-lookup/*.txt);
	do
		sample_name=$(echo ${x##*/} | cut -d'_' -f1)
		tumour_number=$(echo ${x##*/} | cut -d'T' -f1)
		#### check if matched host exists
		if [ -f ${VARIANTS}/host/substitutions/${tumour_number}H* ]
		then
			matched_host=$(ls ${VARIANTS}/host/substitutions/${tumour_number}H*)
			echo -e "\nFiltering $sample_name against ${matched_host} \n" 
			python ${SCRIPTS}/lookuptable.py ${VARIANTS}/tumour/substitutions/pre-lookup/${sample_name}.txt ${matched_host} | sort -k 2,2n > ${VARIANTS}/tumour/substitutions/post-lookup/${sample_name}.txt
			elif [ ! -f ${VARIANTS}/host/substitutions/${tumour_number}H* ];
			then
			 echo -e "\nFiltering unmatched tumour sample ${sample_name} \n"  
			 echo "$sample_name" >> ${SCRIPTS}/unmatched_tumours.txt
			fi
		done
}
