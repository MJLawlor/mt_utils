#!/bin/bash

# MITOTYPUS DATA ANALYSIS PIPELINE

# INPUT
# $1: path to folder containing PLATYPUS output

if [ "$#" -ne 4 ]; then
	echo "| filter_GenotypedVariants.sh: Creates variant lists and VAF plots from Platypus output"
	echo "|"
	echo "| Usage: [filter_GenotypedVariants.sh /path/to/GenotypedVariants_final.vcf] [number of samples included in analysis] [list of unmatched tumours] [list of samples to be discarded]"
	exit	
fi

#### define global variables
PLATYPUSOUT="$1"
PIPELINE="${1}/pipeline"
VARIANTS="${1}/pipeline/variant_lists"
SAMPLES="$2"
SCRIPTS="/lustre/scratch112/sanger/casm/cgp/ml28/MT_2.0/scripts"
LOGS="${1}/pipeline/logs"
DATE=$(date +%y-%m-%d)
REFERENCE="/nfs/dog_n_devil/dog/mitochondria/Canis_familiaris.CanFam3.1.75.dna_sm.chromosome.MT.fa"
mkdir -p ${PIPELINE} ${VARIANTS} ${VARIANTS}/tumour ${VARIANTS}/host ${VARIANTS}/tumour/substitutions ${VARIANTS}/host/substitutions ${VARIANTS}/host/substitutions/pre-filtering ${VARIANTS}/host/substitutions/post-filtering ${VARIANTS}/tumour/sub
stitutions/pre-filtering ${VARIANTS}/host/indels ${VARIANTS}/tumour/indels ${VARIANTS}/tumour/substitutions/pre-lookup ${VARIANTS}/tumour/substitutions/post-lookup ${PIPELINE}/discarded_samples ${PIPELINE}/vaf ${PIPELINE}/trees ${PIPELINE}/trees/
tab ${PIPELINE}/trees/vcf ${PIPELINE}/trees/fasta ${PIPELINE}/trees/fasta/tumour ${LOGS} 
#### create array from list of tumour samples with host contamination
#IFS=$'\n' read -d '' -r -a tum_hostcontam < ${3}
#### create array from list of unmatched tumour samples
IFS=$'\n' read -d '' -r -a unmatched_tumour < ${3}
#### create array from list of samples to be discarded
IFS=$'\n' read -d '' -r -a samples_to_discard < ${4}
#### define function to check if a sample is in an array
function checkArray {
  local i
  for i in "${@:2}"; do [[ "$i" == "$1" ]] && return 0; done
  return 1
}

#### define array of clade defining variants
clade_defining=(229 295 395 530 885 955 1351 1481 2322 2683 2962 3196 3628 3641 3897 3936 4940 5367 5444 5700 6065 6068 6257 6302 6401 6494 6554 6744 7324 7593 8242 8281 8368 8408 8565 8807 8961 9372 9790 9896 10304 10313 10467 10497 10611 10618 
10992 11657 11813 11945 11975 12290 12330 12336 12345 12505 12858 13299 13620 13971 14554 14977 15214 15613 15620 15627 15814 16025)

#### pipeline functions
function makeUnfilteredVariantLists {
	#### extract variants lists from Platypus
	#### sort variant lists into indels and substitutions
	#### discard any variants called in the simple repeats region 16129bp-16430bp across all samples
	echo -e "\n(1) EXTRACTING VARIANT LISTS FROM GenotypedVariants_final.vcf\n"
	input=${PLATYPUSOUT}/GenotypedVariants_final.vcf
	loopvar=$(($SAMPLES + 9))

	egrep -v '#' "$input" > ${PLATYPUSOUT}/data.tmp.vcf
	egrep -v '##' "$input" | head -1  > ${PLATYPUSOUT}/sample_names.tmp.vcf
	for i in $(seq 10 "$loopvar");
	do 
		name=$(awk -v y="$i" '{print $y}' ${PLATYPUSOUT}/sample_names.tmp.vcf | sed 's/^00//g' | sed 's/^0//g' | sed 's/^000//g')
		echo $name
		#name=$(awk -v y="$i" '{print $y}' ${PLATYPUSOUT}/sample_names.tmp.vcf | cut -d'.' -f1)
		echo -e "\nGenerating variant list for $name\n"	
		#### also excluding any variants in the simple repeats region
		awk -v z="$i" '{print $1,$2,$4,$5,$6,$7,$z}' ${PLATYPUSOUT}/data.tmp.vcf |sed -e 's/:/ /g'| awk '{if($12!="0"){$13=$12 / $11 ;print $1,$2,$3">"$4,$5,$6,$13}}' | awk '{if(length($3) == 3){print $0}}' | awk '{if($2<16129 || $2>16430
){print $0}}'| awk '{if($2<15524 || $2>15535){print $0}}'| awk '{if($2<15926 || $2>15938){print $0}}' | awk '{if($2<16662 || $2>16674){print $0}}' | awk '{if($2<1 || $2>75){print $0}}'| awk '{if($2<16652 || $2>16727){print $0}}' > ${VARIANTS}/${n
ame}_substitutions.txt
		awk -v z="$i" '{print $1,$2,$4,$5,$6,$7,$z}' ${PLATYPUSOUT}/data.tmp.vcf |sed -e 's/:/ /g'| awk '{if($12!="0"){$13=$12 / $11 ;print $1,$2,$3">"$4,$5,$6,$13}}' | awk '{if(length($3) > 3){print $0}}' | awk '{if($2<16129 || $2>16430)
{print $0}}' > ${VARIANTS}/${name}_indels.txt
	done
	## sort variant lists into folders
	mv ${VARIANTS}/*H*_indels.txt ${VARIANTS}/host/indels
	mv ${VARIANTS}/*H*_substitutions.txt ${VARIANTS}/host/substitutions/pre-filtering
	mv ${VARIANTS}/*T*_indels.txt ${VARIANTS}/tumour/indels
	mv ${VARIANTS}/*T*_substitutions.txt ${VARIANTS}/tumour/substitutions/pre-filtering
	###clean up
	rm ${PLATYPUSOUT}/*tmp*
}

function fixSampleswap {
	mv ${VARIANTS}/host/substitutions/pre-filtering/1088H-Dog_substitutions.txt ${VARIANTS}/host/substitutions/pre-filtering/1089H-Dog_substitutions.tmp.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/1089H-Dog_substitutions.txt ${VARIANTS}/host/substitutions/pre-filtering/1088H-Dog_substitutions.tmp.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/1089H-Dog_substitutions.tmp.txt ${VARIANTS}/host/substitutions/pre-filtering/1089H-Dog_substitutions.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/1088H-Dog_substitutions.tmp.txt ${VARIANTS}/host/substitutions/pre-filtering/1088H-Dog_substitutions.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/622H-Dog.MT.reheader_substitutions.txt ${VARIANTS}/host/substitutions/pre-filtering/623H-Dog.MT.reheader_substitutions.tmp.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/623H-Dog.MT.reheader_substitutions.txt ${VARIANTS}/host/substitutions/pre-filtering/622H-Dog.MT.reheader_substitutions.tmp.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/622H-Dog.MT.reheader_substitutions.tmp.txt ${VARIANTS}/host/substitutions/pre-filtering/622H-Dog.MT.reheader_substitutions.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/623H-Dog.MT.reheader_substitutions.tmp.txt ${VARIANTS}/host/substitutions/pre-filtering/623H-Dog.MT.reheader_substitutions.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/637H-Dog.reheader_substitutions.txt ${VARIANTS}/host/substitutions/pre-filtering/639H-Dog.reheader_substitutions.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/535H-Dog.reheader_substitutions.txt ${VARIANTS}/host/substitutions/pre-filtering/5035H-Dog.reheader_substitutions.txt 
	mv ${VARIANTS}/host/substitutions/pre-filtering/535H-Dog.MT.reheader_substitutions.txt ${VARIANTS}/host/substitutions/pre-filtering/5035H-Dog.MT.reheader_substitutions.txt 
}

function lookupStep {
	#### compare variants called in matched host tumour samples by PM value and discard host contamination
	#### exclude anything in hosts with PM<0.1
	#### exclude anything in filtered tumour with PM<0.1
	echo -e "\n(3) FILTERING VARIANT LISTS. EXECUTING LOOKUP STEP\n"
	for k in $(ls ${VARIANTS}/host/substitutions/pre-filtering/*.txt);
	do
		sample_name=$(echo ${k##*/} | cut -d'_' -f1)
		echo -e "\nFiltering $sample_name variants list\n" 
		host_number=$(echo ${k##*/} | cut -d'H' -f1)
		awk '{if($6>=0.5){print $0}}'  ${k} > ${VARIANTS}/host/substitutions/post-filtering/${sample_name}.txt
	done
	for j in $(ls ${VARIANTS}/tumour/substitutions/pre-filtering/*.txt);
	do
		sample_name=$(echo ${j##*/} | cut -d'_' -f1)
		tumour_number=$(echo ${j##*/} | cut -d'T' -f1)
		awk '{if($6>=0.2){print $0}}'  ${j} > ${VARIANTS}/tumour/substitutions/pre-lookup/${sample_name}.txt
		#### check if matched host filtered variants list exists
		if [[ $sample_name == *MT* ]]
		then
			if [ -f ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H*MT* ]
			then
			matched_host=$(ls ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H*MT*)
			echo -e "\nFiltering $sample_name against ${matched_host} \n" 
			python ${SCRIPTS}/host_tumour_check.py ${VARIANTS}/tumour/substitutions/pre-lookup/${sample_name}.txt ${matched_host} | sort -k 2,2n > ${VARIANTS}/tumour/substitutions/post-lookup/${sample_name}.txt
			elif [ -f ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H* ]
			then
			matched_host=$(ls ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H*)			
			echo -e "\nFiltering $sample_name against ${matched_host} \n" 	
			python ${SCRIPTS}/host_tumour_check.py ${VARIANTS}/tumour/substitutions/pre-lookup/${sample_name}.txt ${matched_host} | sort -k 2,2n > ${VARIANTS}/tumour/substitutions/post-lookup/${sample_name}.txt
			elif [ ! -f ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H* ];
			then
			echo -e "\nFiltering unmatched tumour sample ${sample_name} \n"  
			#echo "$sample_name" >> ${SCRIPTS}/unmatched_tumours.txt
			awk '{if($6>=0.5){print $0}}' ${VARIANTS}/tumour/substitutions/pre-lookup/${sample_name}.txt > ${VARIANTS}/tumour/substitutions/post-lookup/${sample_name}.txt
			fi
			fi
		
			if [[ $sample_name != *MT* ]]
                	then
                        if [ -f ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H*MT* ]
                        then
                        matched_host=$(ls ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H*MT*)
                        echo -e "\nFiltering $sample_name against ${matched_host} \n" 
                        python ${SCRIPTS}/host_tumour_check.py ${VARIANTS}/tumour/substitutions/pre-lookup/${sample_name}.txt ${matched_host} | sort -k 2,2n > ${VARIANTS}/tumour/substitutions/post-lookup/${sample_name}.txt
                        elif [ -f ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H* ]
                        then
                        matched_host=$(ls ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H*)
                        echo -e "\nFiltering $sample_name against ${matched_host} \n"    
                        python ${SCRIPTS}/host_tumour_check.py ${VARIANTS}/tumour/substitutions/pre-lookup/${sample_name}.txt ${matched_host} | sort -k 2,2n > ${VARIANTS}/tumour/substitutions/post-lookup/${sample_name}.txt
                        elif [ ! -f ${VARIANTS}/host/substitutions/post-filtering/${tumour_number}H* ];
                        then
                        echo -e "\nFiltering unmatched tumour sample ${sample_name} \n"  
                	#echo "$sample_name" >> ${SCRIPTS}/unmatched_tumours.txt
                	awk '{if($6>=0.5){print $0}}' ${VARIANTS}/tumour/substitutions/pre-lookup/${sample_name}.txt > ${VARIANTS}/tumour/substitutions/post-lookup/${sample_name}.txt
                	fi
                	fi
		done
}
function rescueVariants {
	echo -e "\n(3) RECOVERING DISCARDED VARIANTS IN TUMOURS WITH TUMOUR CONTAMINATED HOST SAMPLES\n"
#	mkdir -p ${VARIANTS}/tumour/substitutions/post-filtering/variant_lists_before_rescue
	#### rescue real variants discarded from tumours with tumour contaminated hosts
		for i in ${tum_hostcontam[@]};
		do
			echo -e "\nRescuing variants in $i\n"
			python ${SCRIPTS}/rescueVariants.py ${VARIANTS}/tumour/substitutions/pre-filtering/${i}* ${VARIANTS}/tumour/substitutions/post-filtering/${i}* ${VARIANTS}/tumour/substitutions/post-filtering/${i}.rescue.tmp.txt 
			mv ${VARIANTS}/tumour/substitutions/post-filtering/${i}.txt ${VARIANTS}/tumour/substitutions/post-filtering/variant_lists_before_rescue
			sort -k 2,2n ${VARIANTS}/tumour/substitutions/post-filtering/${i}.rescue.tmp.txt | sed 's/\t/ /g' > ${VARIANTS}/tumour/substitutions/post-filtering/${i}.txt 
			rm -r ${VARIANTS}/tumour/substitutions/post-filtering/${i}.rescue.tmp.txt 
		done 
}

function addHeader {
	echo "CHROM POS REF>ALT QUAL FILTER PM" > ${PLATYPUSOUT}/header.tmp.txt
	echo -e "\n(4) ADDING HEADERS TO VARIANT LISTS\n"
	mkdir -p ${VARIANTS}/${1}/substitutions/${2}/headers
	for j in $(ls ${VARIANTS}/${1}/substitutions/${2}/*Dog*);
	do
		path=$(dirname "${j}")
		sample_name=$(echo ${j##*/} | sed 's/.txt//g'|cut -d'_' -f1)
		echo -e "\nAdding header to ${sample_name}\n"
		cat ${PLATYPUSOUT}/header.tmp.txt $j > ${path}/headers/${sample_name}.header.txt
	done
	rm ${PLATYPUSOUT}/*tmp*
}


function generateVAFplots {
	echo -e "\n(5) CREATING VAF PLOTS FOR TUMOUR AND HOST SAMPLES\n"
	for j in $(ls ${VARIANTS}/${1}/substitutions/pre-filtering/headers/*Dog*txt);
	do	
	text_file=$(basename "${j}")
	tumour_name=$(echo $text_file | cut -d'-' -f1)
	sample_name=$(echo $text_file | sed 's/.txt//g')
	sample=$(echo $sample_name | sed 's/\.header//g')
	path=$(dirname "${j}")
	checkArray "${tumour_name}-Dog" "${unmatched_tumour[@]}"
	check=$(echo $?)
	if [ $check == "1" ] && [[ $tumour_name == *"T"* ]];
	then
	echo -e "\nGenerating VAF plot for matched tumour${tumour_name}\n"
	Rscript $SCRIPTS/vaf.r $j ${VARIANTS}/${1}/substitutions/post-lookup/headers/$text_file ${PIPELINE}/vaf/${sample_name}.pdf $sample
	Rscript $SCRIPTS/vaf.r $j ${VARIANTS}/${1}/substitutions/post-lookup/headers/$text_file ${PIPELINE}/vaf/${sample_name}_labelled.pdf $sample label
	elif [[ $tumour_name == *"H"* ]]
	then
	echo -e "\nGenerating VAF plot for host${tumour_name}\n"
	Rscript $SCRIPTS/vaf.r $j ${VARIANTS}/${1}/substitutions/post-filtering/headers/$text_file ${PIPELINE}/vaf/${sample_name}.pdf $sample 
	Rscript $SCRIPTS/vaf.r $j ${VARIANTS}/${1}/substitutions/post-filtering/headers/$text_file ${PIPELINE}/vaf/${sample_name}_labelled.pdf $sample label
	elif [ $check == "0" ]
	then
	echo -e "\nGenerating VAF plot for unmatched tumour${tumour_name}\n"	
	Rscript $SCRIPTS/vaf.r ${j} ${VARIANTS}/${1}/substitutions/post-lookup/headers/$text_file ${PIPELINE}/vaf/${sample_name}.pdf $sample 
	Rscript $SCRIPTS/vaf.r ${j} ${VARIANTS}/${1}/substitutions/post-lookup/headers/$text_file ${PIPELINE}/vaf/${sample_name}_labelled.pdf $sample label 
	fi
	done

	#gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="${PIPELINE}/vaf/${DATE}_tumours.pdf" ${PIPELINE}/vaf/*T*-Dog.pdf
	#gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="${PIPELINE}/vaf/${DATE}_hosts.pdf" ${PIPELINE}/vaf/*H*-Dog.pdf
	#gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="${PIPELINE}/vaf/${DATE}_tumours_labelled.pdf" ${PIPELINE}/vaf/*T*-Dog_labelled.pdf
	#gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="${PIPELINE}/vaf/${DATE}_hosts_labelled.pdf" ${PIPELINE}/vaf/*H*-Dog_labelled.pdf

###for low resolution pdfs use
#gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=
}

function assignClades {
	echo -e "\n(1) ASSIGNING TUMOUR SAMPLES TO CLADES \n"
	for j in $(ls ${VARIANTS}/tumour/substitutions/post-lookup/*.txt);
	do
	python ${SCRIPTS}/clade.py $j >> ${SCRIPTS}/${DATE}_clades.txt
	done
}

function discardSamples {
	echo -e "\n(1) DISCARDING SAMPLES \n"
	for i in $(ls ${VARIANTS}/${1}/substitutions/pre-filtering/*Dog*);
	do 
		sample_name=$(echo ${i##*/} | cut -d'_' -f1)
		checkArray "${sample_name}" "${samples_to_discard[@]}"
       	 	check=$(echo $?)
	        if [ $check == "0" ];
		then
			echo -e "\nDiscarding sample ${sample_name}\n" 
			rm ${i}
			###rm ${VARIANTS}/${1}/substitutions/pre-lookup/${sample_name}.txt 
			###mv ${VARIANTS}/${1}/substitutions/post-lookup/${sample_name}.txt ${PIPELINE}/discarded_samples/${sample_name}.post-filtering.txt
		####elif [ $check == "0" ] && [[ $sample_name == *"H"* ]];
		###then
		###	echo -e "\nDiscarding sample ${sample_name}\n" 
		###	rm ${i}
		###	mv ${VARIANTS}/${1}/substitutions/post-filtering/${sample_name}.txt ${PIPELINE}/discarded_samples/${sample_name}.txt
		fi
	done
}
function makeTrees {
 	#converts tab delimited file into a vcf and then into a fasta file
        #input tab has to be in the format CHROM        POS     RefAlt
	echo -e "\n(1) CREATING PHYLOGENETIC TREE FROM FINALIZED VARIANT LISTS \n"
        for i in ${VARIANTS}/host/substitutions/post-filtering/*.txt;
#        for i in ${VARIANTS}/tumour/substitutions/post-lookup/*.txt;
	do
			sample_name=$(echo ${i##*/} | sed 's/.txt//g')
			echo "Making tab file for $sample_name"
			awk '{print $1"\t"$2"\t"$3}' ${i} | sed 's/>/ /g'| sort -k 2,2n > ${PIPELINE}/trees/tab/${sample_name}.tab
#                        #tail --lines=+2 ${trees}/${w} | awk '{print $1"\t"$2"\t"$3$4}' > ${trees}/${y}.tab
                        echo "Success!"
                        echo "Making vcf from ${sample_name} tab file"
                        rm -f ${PIPELINE}/trees/vcf/${sample_name}.vcf
                        #using most recent verision of tab-to-vcf at /software/CGP
                       	#cat ${PIPELINE}/trees/tab/${sample_name}.tab | tab-to-vcf -i ${PIPELINE}/trees/tab/${sample_name} -r $REFERENCE > ${PIPELINE}/trees/vcf/${sample_name}.vcf
			echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	${sample_name}" > ${PIPELINE}/trees/vcf/header.tmp
			awk '{print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\t.\tGT\t0/1"}' ${PIPELINE}/trees/tab/${sample_name}.tab > ${PIPELINE}/trees/vcf/${sample_name}.tmp
			cat ${PIPELINE}/trees/vcf/header.txt ${PIPELINE}/trees/vcf/header.tmp ${PIPELINE}/trees/vcf/${sample_name}.tmp > ${PIPELINE}/trees/vcf/${sample_name}.vcf
			echo "Success!"
                        bgzip -f ${PIPELINE}/trees/vcf/${sample_name}.vcf
                       	/software/CGP/bin/tabix -p vcf ${PIPELINE}/trees/vcf/${sample_name}.vcf.gz
                        echo "Making fasta file from ${sample_name} vcf"
			cat $REFERENCE | vcf-consensus ${PIPELINE}/trees/vcf/${sample_name}.vcf.gz |sed -e "s/>MT/>"${sample_name}"/g" > ${PIPELINE}/trees/fasta/tumour/${sample_name}.fa
			echo "Success!"
        done
       cat ${PIPELINE}/trees/fasta/tumour/*.fa > ${PIPELINE}/trees/fasta/${DATE}_tumour_samples.fa
       #perl ${SCRIPTS}/Fasta2Phylip.pl ${PIPELINE}/trees/fasta/${DATE}_tumour_samples.fa ${PIPELINE}/trees/phy/${DATE}_tumour_samples.phy
}
function runPhyml() {
	bsub -R"select[mem>15900] rusage[mem=15900]" -M15900 -o ${LOGS}/${DATE}.%J.stdout -e ${LOGS}/${DATE}.%J.stderr phyml -i ${1} -b -4 -m GTR -t e -v e -a e -s BEST -o tlr --run_id=${1} --quiet
}

# List of available functions
#makeUnfilteredVariantLists
#fixSampleswap
#discardSamples tumour
#discardSamples host
#lookupStep
#rescueVariants
#addHeader host pre-filtering
#addHeader host post-filtering
#addHeader tumour pre-filtering
#addHeader tumour post-lookup
#generateVAFplots host
#generateVAFplots tumour
#assignClades
#makeTrees
#runPhyml ${PIPELINE}/trees/phy/${DATE}_tumour_samples.fa
