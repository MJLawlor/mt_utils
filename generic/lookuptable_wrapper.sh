#!/bin/bash

# mitotypus data post-processing pipeline: lookup step wrapper script
# author:       Máire Ní Leathlobhair (ml677@cam.ac.uk)
# date:         May 2016

# input
# $1: absolute path to folder containing Mitotypus output
# $2: total number of samples (including host and tumour)
# $3: absolute path to scripts folder
# $4: absolute path to MT reference genome
# $5: VAF value used in tumour vs. matched normal lookup step
# $6: path to raxml

if [ "$#" -ne 6 ]; then
	echo "|"
	echo "| lookuptable_wrapper.sh: Creates variant lists and VAF plots from Mitotypus output"
	echo "|"
	echo "| Usage: lookuptable_wrapper.sh /path/to/GenotypedVariants_final.vcf [number of samples to be analysed] /path/to/scripts/folder /path/to/reference/MT/genome VAF/cutoff /path/to/raxml/install"
	echo "|"
	exit
fi

PLATYPUSOUT="$1"
PIPELINE="${1}/pipeline"
VARIANTS="${1}/pipeline/variant_lists"
SAMPLES="$2"
SCRIPTS="$3"
LOGS="${1}/pipeline/logs"
DATE=$(date +%y-%m-%d)
REFERENCE="$4"
CUTOFF="$5"
RAXML="$6"
mkdir -p ${PIPELINE} ${LOGS} ${VARIANTS} ${VARIANTS}/tumour ${VARIANTS}/host ${VARIANTS}/tumour/substitutions ${VARIANTS}/host/substitutions ${VARIANTS}/host/indels ${VARIANTS}/tumour/indels ${VARIANTS}/tumour/substitutions/pre-lookup ${VARIANTS}/tumour/substitutions/post-lookup ${PIPELINE}/vaf 

function checkArray {
  local i
  for i in "${@:2}"; do [[ "$i" == "$1" ]] && return 0; done
  return 1
}

function rawVariantlists { 
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
		awk -v z="$i" '{print $1,$2,$4,$5,$6,$7,$z}' ${PLATYPUSOUT}/data.tmp.vcf |sed -e 's/:/ /g'| awk '{if($12!="0"){$13=$12 / $11 ;print $1,$2,$3">"$4,$13}}' | awk '{if(length($3) == 3){print $0}}' > ${VARIANTS}/${name}_substitutions.txt
		awk -v z="$i" '{print $1,$2,$4,$5,$6,$7,$z}' ${PLATYPUSOUT}/data.tmp.vcf |sed -e 's/:/ /g'| awk '{if($12!="0"){$13=$12 / $11 ;print $1,$2,$3">"$4,$13}}' | awk '{if(length($3) > 3){print $0}}' > ${VARIANTS}/${name}_indels.txt
	done
	mv ${VARIANTS}/*H*_indels.txt ${VARIANTS}/host/indels
	mv ${VARIANTS}/*H*_substitutions.txt ${VARIANTS}/host/substitutions/
	mv ${VARIANTS}/*T*_indels.txt ${VARIANTS}/tumour/indels
	mv ${VARIANTS}/*T*_substitutions.txt ${VARIANTS}/tumour/substitutions/pre-lookup
	rm ${PLATYPUSOUT}/*tmp*
}

function lookupStep {
# compare variants called in matched host tumour samples by PM value and discard host contamination
 echo -e "\n(2) FILTERING VARIANT LISTS. EXECUTING LOOKUP STEP\n"
	rm -f ${SCRIPTS}/unmatched_tumours.txt # prevent continuously writing to this text file when re-running wrapper script
	for x in $(ls ${VARIANTS}/tumour/substitutions/pre-lookup/*.txt);
	do
		sample_name=$(echo ${x##*/} | cut -d'_' -f1)
		tumour_number=$(echo ${x##*/} | cut -d'T' -f1)
		# check if matched host exists
		if [ -f ${VARIANTS}/host/substitutions/${tumour_number}H* ]
		then
			matched_host=$(ls ${VARIANTS}/host/substitutions/${tumour_number}H*)
			echo -e "\nFiltering $sample_name against ${matched_host} \n" 
			python ${SCRIPTS}/lookuptable.py ${VARIANTS}/tumour/substitutions/pre-lookup/${sample_name}_substitutions.txt ${matched_host} | sort -k 2,2n > ${VARIANTS}/tumour/substitutions/post-lookup/${sample_name}.txt
			elif [ ! -f ${VARIANTS}/host/substitutions/${tumour_number}H* ];
			then
			 echo -e "\nFiltering unmatched tumour sample ${sample_name} \n"  
			 echo "$sample_name" >> ${SCRIPTS}/unmatched_tumours.txt
			fi
		done
}

function vafPlot {
	echo "CHROM POS REF>ALT VAF" > ${PLATYPUSOUT}/header.tmp.txt
        echo -e "\n(3) ADDING HEADERS TO VARIANT LISTS\n"
        mkdir -p ${VARIANTS}/tumour/substitutions/pre-lookup/headers ${VARIANTS}/tumour/substitutions/post-lookup/headers ${VARIANTS}/host/substitutions/headers
        for j in $(ls ${VARIANTS}/tumour/substitutions/pre-lookup/*.txt ${VARIANTS}/tumour/substitutions/post-lookup/*.txt ${VARIANTS}/host/substitutions/*.txt);
        do
                path=$(dirname "${j}")
                sample_name=$(echo ${j##*/} | sed 's/.txt//g'|cut -d'_' -f1)
                echo -e "\nAdding header to ${sample_name}\n"
                cat ${PLATYPUSOUT}/header.tmp.txt $j > ${path}/headers/${sample_name}.header.txt
        done
        rm ${PLATYPUSOUT}/*tmp*
        echo -e "\n(4) CREATING VAF PLOTS FOR TUMOUR AND HOST SAMPLES\n"
        for x in $(ls ${VARIANTS}/tumour/substitutions/pre-lookup/headers/*txt ${VARIANTS}/host/substitutions/headers/*txt);
        do
         text_file=$(basename "${x}")
         tumour_name=$(echo $text_file | cut -d'-' -f1)
         sample_name=$(echo $text_file | sed 's/.txt//g')
         sample=$(echo $sample_name | sed 's/\.header//g')
         path=$(dirname "${x}")
         IFS=$'\n' read -d '' -r -a unmatched_tumour < ${SCRIPTS}/unmatched_tumours.txt
         checkArray "${tumour_name}-Devil" "${unmatched_tumour[@]}"
         check=$(echo $?)
         if [ $check == "1" ] && [[ $tumour_name == *"T"* ]];
         then
          echo -e "\nGenerating VAF plot for matched tumour ${tumour_name}\n"
          Rscript $SCRIPTS/vaf.R $x ${VARIANTS}/tumour/substitutions/post-lookup/headers/$text_file ${PIPELINE}/vaf/${sample_name}.pdf $CUTOFF $sample
          Rscript $SCRIPTS/vaf.R $x ${VARIANTS}/tumour/substitutions/post-lookup/headers/$text_file ${PIPELINE}/vaf/${sample_name}_labelled.pdf $CUTOFF $sample label
         elif [[ $tumour_name == *"H"* ]]
         then
          echo -e "\nGenerating VAF plot for host ${tumour_name}\n"
          Rscript $SCRIPTS/vaf.R $x ${VARIANTS}/host/substitutions/headers/$text_file ${PIPELINE}/vaf/${sample_name}.pdf $CUTOFF $sample
          Rscript $SCRIPTS/vaf.R $x ${VARIANTS}/host/substitutions/headers/$text_file ${PIPELINE}/vaf/${sample_name}_labelled.pdf $CUTOFF $sample label
         elif [ $check == "0" ]
         then
          echo -e "\nGenerating VAF plot for unmatched tumour ${tumour_name}\n"    
          Rscript $SCRIPTS/vaf.R ${x} ${x} ${PIPELINE}/vaf/${sample_name}.pdf $CUTOFF $sample
          Rscript $SCRIPTS/vaf.R ${x} ${x} ${PIPELINE}/vaf/${sample_name}_labelled.pdf $CUTOFF $sample label
         fi
        done
}

function phylogeny {
 	# converts tab delimited variant list into vcf format and then into fasta format
        # input tab has to be in the format CHROM        POS     RefAlt
        mkdir ${PIPELINE}/trees ${PIPELINE}/trees/vcf ${PIPELINE}/trees/phy ${PIPELINE}/trees/tab ${PIPELINE}/trees/fasta ${PIPELINE}/trees/raxml
        # create generic vcf header
        printf "##fileformat=VCFv4.0 \n##fileDate=%s\n" "$DATE" > ${PIPELINE}/trees/vcf/header.txt
	echo -e "\n(5) CREATING FASTAS FROM FINALIZED VARIANT LISTS \n"
        for i in $(ls ${VARIANTS}/host/substitutions/*.txt ${VARIANTS}/tumour/substitutions/post-lookup/*.txt );do
	 sample_name=$(echo ${i##*/} | sed 's/.txt//g')
	 echo -e "\nMaking tab file for ${sample_name}\n"
	 awk '{print $1"\t"$2"\t"$3}' ${i} | sed 's/>/ /g'| sort -k 2,2n > ${PIPELINE}/trees/tab/${sample_name}.tab
         echo "Making vcf from ${sample_name} tab file"
         rm -f ${PIPELINE}/trees/vcf/${sample_name}.vcf
         # using most recent verision of tab-to-vcf at /software/CGP
	 echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	${sample_name}" > ${PIPELINE}/trees/vcf/header.tmp
	 awk '{print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\t.\tGT\t0/1"}' ${PIPELINE}/trees/tab/${sample_name}.tab > ${PIPELINE}/trees/vcf/${sample_name}.tmp
	 cat ${PIPELINE}/trees/vcf/header.txt ${PIPELINE}/trees/vcf/header.tmp ${PIPELINE}/trees/vcf/${sample_name}.tmp > ${PIPELINE}/trees/vcf/${sample_name}.vcf
         bgzip -f ${PIPELINE}/trees/vcf/${sample_name}.vcf
         /software/CGP/bin/tabix -p vcf ${PIPELINE}/trees/vcf/${sample_name}.vcf.gz
         echo -e "\nMaking fasta file from ${sample_name}.vcf\n"
	 cat $REFERENCE | vcf-consensus ${PIPELINE}/trees/vcf/${sample_name}.vcf.gz |sed -e "s/>MT/>"${sample_name}"/g" > ${PIPELINE}/trees/fasta/${sample_name}.fa
        done
       cat ${PIPELINE}/trees/fasta/*.fa > ${PIPELINE}/trees/raxml/${DATE}_all_samples.fa
       perl ${SCRIPTS}/Fasta2Phylip.pl ${PIPELINE}/trees/raxml/${DATE}_all_samples.fa ${PIPELINE}/trees/raxml/${DATE}_all_samples.phy
}
function runRaxml() {
	# run RAxML using the GTR+GAMMA+P-Invar model and default rapid hill climbing algorithm
	echo -e "\n(6) RUNNING PHYLOGENETIC ANALYSIS USING RAxML\n"
	bsub -R"select[mem>15900] rusage[mem=15900]" -M15900 -o ${LOGS}/${DATE}_RAxML.%J.stdout -e ${LOGS}/${DATE}_RAxML.%J.stderr ${RAXML}raxmlHPC -m GTRGAMMAI -s ${PIPELINE}/trees/raxml/${DATE}_all_samples.phy -n ${DATE} -p $RANDOM
}

rawVariantlists
lookupStep
vafPlot
phylogeny
runRaxml
