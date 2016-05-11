#!/bin/bash
# Last update: 09/12/14 by Maire

# globs
today=`date +%Y-%m-%d`
runHOME="$1" # set home directory
HOST="${runHOME}/Host" # create sub-folder for host samples
TUMOUR="${runHOME}/Tumour" # create sub-folder for tumour samples
PINDELCGP="$2" # path for pindel output
CAVEMAN="$3" # path for caveman output
PIPELINERESULTS="/nfs/dog_n_devil/dog/mitochondria/mt_pipeline"
RESULTSMUTS="${PIPELINERESULTS}/caveman_output/MUTS"
RESULTSSNPS="${PIPELINERESULTS}/caveman_output/SNPS"
RESULTSPINDEL="${PIPELINERESULTS}/pindel_output"
RESULTSMPILE="${PIPELINERESULTS}/mpileup_output"
SCRIPTS="${PIPELINERESULTS}/analysis/scripts"
POLYC="${runHOME}/polyc_graphs"
BED="${runHOME}/bed"

# make sure all of these folders exist
	mkdir -p "$TUMOUR" "$HOST"

#  functions
function canapps_download()
{
	bams_to_download=($(awk '{print $2}' sample_names_list.txt | grep ${1}))
	
	for a in "${bams_to_download[@]}";
		do
		# if bam file does not exist then download it else print message to user
		if [ ! -f ${2}/${a}/${a}.bam ]; then
			download_bam=$(pipelineResult -db live -p 1115 -s $a -t 7)
			echo "Downloading $a from $download_bam"
			# check file size before and after
			ls -lh $download_bam | awk '{print $5}'
			if [ ! -d ${2}/${a}/ ];then
				echo "Making folder for sample ${a}."
				mkdir -p ${2}/${a}/
			else
				echo "Folder for sample ${a}/ already exists."
			fi
			cp $download_bam ${2}/${a}/"$a".bam
			ls -lh ${2}/${a}/"$a".bam | awk '{print $5}'
			# download bam file index
			download_index=$(pipelineResult -db live -p 1115 -s $a -t 16)
			echo "Downloading $a from $download_index"
			ls -lh $download_index | awk '{print $5}'
			cp $download_index ${2}/${a}/"$a".bam.bai
			ls -lh ${2}/${a}/"$a".bam.bai | awk '{print $5}'
		else 
			echo "Bam file ${a}.bam  already exists."
		fi
		done
}

function check_file_size()
{
  # moves bam files to failed folder if they fail the file size check
  echo "Checking bam file sizes"
	for b in $(ls -d "$TUMOUR"/* "$HOST"/*);do
		b_samplename=$(echo ${b##*/} | cut -d'.' -f1)
		if [ -f ${b}/${b_samplename}.bam ];then
			limit=70000000 
			size=$(ls -l ${b}/${b_samplename}.bam | awk '{print $5}')
	       		if [ "$size" -lt "$limit" ];then
                		echo "${b_samplename} failed file size check. $size < $limit"
				mv ${b}/${b_samplename}.bam "$runHOME"/failed
			else
				echo "${b_samplename} is fine. $size > $limit"
			fi
        	else
			echo "${b_samplename}.bam does not exist"
		fi
 done
}

function index_sort()
{
	# index and sort bam file. remove duplicates.
    for j in $(ls -d "$@");do
      ij=$(echo ${j##*/} | cut -d'.' -f1) #name of sample
		  # make sure bam file already exists
		  if [ -f $j/${ij}.bam ];then
			# if folder with sample name does not already exist, make it
			  if [ ! -f $j/${ij}.bam.bai ];then
				  samtools index $j/$ij.bam
			  else 
				  echo "${ij}.bam.bai already exists."
			  fi
	 		# sort bam
			  if [ ! -f $j/$ij.sorted.bam ];then
				  samtools sort $j/${ij}.bam $j/$ij.sorted
			  else
				  echo "${ij}.sorted.bam already exists. "
        fi
			# remove duplicates
			  if [ ! -f $j/$ij.sorted.rmdup.bam ];then
				  samtools rmdup $j/$ij.sorted.bam $j/$ij.sorted.rmdup.bam
			  else
				  echo "$ij.sorted.rmdup.bam already exists."				
			  fi
			# index again
			 if [ ! -f $j/${ij}.sorted.rmdup.bam.bai ];then
        samtools index $j/$ij.sorted.rmdup.bam
       else
				echo "$ij.sorted.rmdup.bam.bai already exists."
			fi
		  else
			 echo "${ij}.bam does not exist"
		 fi
	done
}
function pindel_bams()
{		
	# extract only MT reads from all host and tumour bam files
	for k in $(ls -d "$@");do
		jk=$(echo ${k##*/} | cut -d'.' -f1) # strip path
		if [ -f $k/${jk}.sorted.rmdup.bam ];then
			cd $k
			if [ ! -f $k/$jk.mt_only.sam ];then
				echo "Making MT only bam file for $k"
				samtools view -h $jk.sorted.rmdup.bam MT > $jk.mt_only.sam
			else
				echo "$jk.mt_only.sam already exists."
			fi
			if [ ! -f $k/$jk.mt_only.bam ];then
				samtools view -bS $jk.mt_only.sam -o $jk.mt_only.bam
			else
				echo "$jk.mt_only.bam already exists."
			fi
			if [ ! -f $k/$jk.mt_only.bam.bai ];then
				samtools index $jk.mt_only.bam
			else
				echo "$jk.mt_only.bam.bai already exists."
			fi
		else
			echo "${jk}.sorted.rmdup.bam does not exist"

		fi

	done	 
}
function check_coverage(){
#	for l in $(ls -d "$TUMOUR"/* "$HOST"/*);do
		l_samplename=$(echo ${l##*/} | cut -d'.' -f1)
		echo "${l_samplename}" > ${runHOME}/${today}_sample_coverage.txt
		coverageBed -abam $l/"$l_samplename".sorted.rmdup.bam -b ${runHOME}/bed/MT.bed >> ${runHOME}/${today}_sample_coverage.txt
	done
}

#function check_depth()
#{
	#creates text file showing average depth of each muts.vcf
	#rm -f $runHOME/*depth_flag.txt
	#for l in $(ls -d "$TUMOUR"/* "$HOST"/*);do
		#printf "${l##*/}" >> $runHOME/"$today"_depth_flag.txt
		#egrep -v "^#" $l/results/MT/1_16727.muts.vcf |awk '{split($8,a,";"); print a[1]}'|sed "s/[^0-9]//g"|awk  '{sum+=$1} END {print " "sum/FNR}' >> $runHOME/"$today"_depth_flag.txt
		#done
#}

function list_unmatched()
{
	for q in $(ls -d "$TUMOUR"/*); do 
	echo ${q##*/} |cut -d'T' -f1 >> tumour_samples_tmp.txt
	done

	for r in $(ls -d "$HOST"/*); do
        echo ${r##*/} |cut -d'H' -f1 >> host_samples_tmp.txt
        done

 	awk 'FNR==NR {a[$1]; next} !($1 in a) {print $1"T-Dog"}' host_samples_tmp.txt tumour_samples_tmp.txt > "$today"_unmatched_tumour.txt
	awk 'FNR==NR {a[$1]; next} ($1 in a) {print $1"T-Dog"}' host_samples_tmp.txt tumour_samples_tmp.txt > "$today"_matched_tumour.txt
	awk 'FNR==NR {a[$1]; next} !($1 in a) {print $1"H-Dog"}' tumour_samples_tmp.txt host_samples_tmp.txt > "$today"_unmatched_host.txt
	rm -f *tmp.txt

}
function run_pindelcgp()
{
        #run pindelcgp on all host and tumour samples
        for n in $(ls -d "$@");do
		ns=${n##*/} #sample name without path
		if [ -f $n/${ns}.mt_only.bam ];then 
			if [ ! -f ${PINDELCGP}/${ns}/${ns}_vs_DDv3_1is.flagged.vcf ];then	
				cd $n
				echo "Running pindelcgp on ${ns}"
				perl /lustre/scratch104/sanger/ml28/bin/pindel.pl -o /lustre/scratch104/sanger/ml28/pindelcgp/${ns}/ -r "$runHOME"/genome.fa -t $n/${ns}.mt_only.bam -n "$runHOME"/DDv3_1is.mt_only.bam -s "$runHOME"/realsimple
repeats.gff.gz.tbi -f "$runHOME"/ctvtRules.lst -g "$runHOME"/coding_footprint.gff.gz.tbi -u "$runHOME"/unmatched_normal.gff.gz.tbi
        		else
				echo "${ns}_pindel.bed already exists."
			fi
		else
			echo "${ns} does not exist"
		fi
		done

}
function make_pindel_bedfile()
{
	#make germline indel bed file from pindel output. This is used as input for caveman flagging.
	for p in $(ls -d "$PINDELCGP"/*-Dog); do
		samplename=$(echo ${p##*/})
		if [ ! -f ${BED}/${samplename}_indel.bed ];then
			echo "Making germline indel bed file for sample ${samplename}"
	        	if [ -f ${p}/${samplename}_vs_DDv3_1is.flagged.vcf.gz ];then 
				gunzip ${p}/${samplename}_vs_DDv3_1is.flagged.vcf.gz
			fi
			#set coordinates of bed file to be 10bp upstream and downstream of indel offset by REP value
			egrep -v "^#" ${p}/${samplename}_vs_DDv3_1is.flagged.vcf |awk '{split($8,a,";"); a[7]=substr(a[7],5);print "MT""\t"$2-5"\t"$2+5+a[7]}' > ${p}/${samplename}_indel.bed
			cp ${p}/${samplename}_indel.bed ${BED}
		else
			echo "${samplename}.bed already exists."
		fi
	done

}
function run_cavemanwrapper()
{
#run cavemanwrapper script on all bam files which have been properly sorted, indexed etc. Pindel should already have run on these files.
for i in $(ls -d ${1}/807H-Dog);do
	samplename=$(echo ${i##*/})
	if [ -f ${i}/${samplename}.sorted.rmdup.bam ];then
		mkdir -p ${i}/run1 ${i}/run2
		if [ ! -f ${i}/${samplename}_merged.muts.vcf ];then
			echo "Running caveman.pl on ${samplename} using ${samplename}_pindel.bed, ${3}"
			perl /nfs/users/nfs_m/ml28/cgpCaVEManWrapper-1.2.0/bin/caveman.pl -o ${i}/run1 -r "$runHOME"/genome.fa.fai -tb ${2}/${samplename}/${samplename}.sorted.rmdup.bam -nb "$runHOME"/DDv3_1is.sorted.bam -ig "$runHOME"/ignor
e_contigs.tsv -tc "$runHOME"/${3}  -nc "$runHOME"/no_analysis_fake.bed -s DOG -sa CanFam3.1 -b "$runHOME"/bed -in "$runHOME"/bed/${samplename}_indel.bed -u "$runHOME"/empty_vcf -st genomic -l 1 -k 0.0 -c "$runHOME"/flag1.vcf.config.ini
			echo "Finished run."
		else
			echo "Already ran caveman on ${i##*/}"
		fi	
		if [ ! -f ${i}/${i##*/}_merged.muts.vcf ];then
			echo "Running caveman.pl a second time on ${i##*/} using ${k}_pindel.bed, ${3}"
			perl /nfs/users/nfs_m/ml28/cgpCaVEManWrapper-1.2.0/bin/caveman.pl -o ${i}/run2 -r "$runHOME"/genome.fa.fai -tb ${2}/${i##*/}/${i##*/}.sorted.rmdup.bam -nb "$runHOME"/DDv3_1is.sorted.bam -ig "$runHOME"/ignore_contigs.
tsv -tc "$runHOME"/${3}  -nc "$runHOME"/no_analysis_fake.bed -s DOG -sa CanFam3.1 -b "$runHOME"/bed -in "$runHOME"/bed/${i##*/}_indel.bed -u "$runHOME"/empty_vcf -st genomic -l 1 -k 0.0 -c "$runHOME"/flag2.vcf.config.ini
			echo "Finished run."
		else 
			echo "Already ran caveman on ${i##*/}"
		fi
		if [ ! -f ${i}/${i##*/}_merged.muts.vcf ];then
			echo "Merging caveman output for ${i##*/}"
			if [ -f ${i}/run1/${i##*/}_vs_DDv3_1is.flagged.muts.vcf.gz ];then
			gunzip ${i}/run1/${i##*/}_vs_DDv3_1is.flagged.muts.vcf.gz
			fi
			if [ -f ${i}/run2/${i##*/}_vs_DDv3_1is.flagged.muts.vcf.gz ];then
			gunzip ${i}/run2/${i##*/}_vs_DDv3_1is.flagged.muts.vcf.gz
			fi
			#check if host or tumour and rescue variants accordingly
			if [[ ${i##*/} == *H* ]];then
				python "$SCRIPTS"/merge.py ${i}/run1/*.flagged.muts.vcf ${i}/run2/*.flagged.muts.vcf ${i}/${i##*/}_merged.muts.vcf 2683 3028 6882 7014 8281 13708 16663
			elif [[ ${i##*/} == *T* ]];then
				python "$SCRIPTS"/merge.py ${i}/run1/*.flagged.muts.vcf ${i}/run2/*.flagged.muts.vcf ${i}/${i##*/}_merged.muts.vcf 381 1481 1683 2682 6882 8368 8703 9825 9896 14977 15524 15526 16660 16663 16671
			else
				echo "There seems to be a problem. ${i##*/} does not appear to be host or tumour"
			fi
		else
			echo "${i##*/}_merged.muts.vcf already exists."
		fi
else
	echo "${i##*/}.sorted.rmdup.bam does not exist"
        fi
	done
}
function run_mpileup()
{
	#for s in $(ls -d "$TUMOUR"/* "$HOST"/*);do
	for s in $(ls -d "$TUMOUR"/4T-Dog "$HOST"/4H-Dog "$TUMOUR"/2T-Dog "$TUMOUR"/3T-Dog);do
	#run samtools mpileup and exclude indels
		if [ ! -f $s/${s##*/}.raw.bcf ];then
		echo "Running mpileup on ${s##*/}"
		samtools mpileup -f "$runHOME"/genome.fa -l "$runHOME"/bed/MT.bed -g -I ${s}/${s##*/}.sorted.rmdup.bam | bcftools view -bvcg - > ${s}/${s##*/}.raw.bcf
		else
		echo "${s##*/}.raw.bcf already exists."
		fi
		#run samtools mpileup and include indels
		if [ ! -f $s/${s##*/}_plusindels.raw.bcf ];then
                echo "Running mpileup on ${s##*/}"
                samtools mpileup -f "$runHOME"/genome.fa -l "$runHOME"/bed/MT.bed -g ${s}/${s##*/}.sorted.rmdup.bam | bcftools view -bvcg - > ${s}/${s##*/}_plusindels.raw.bcf
                else
                echo "${s##*/}_plusindels.raw.bcf already exists."
                fi
	done
}
function move_caveman_results()
{
        for y in $(ls -d ${1}/*)
                do
	#		if [ -f $y/${y##*/}.raw.bcf ];then
	#        	        cp ${y}/${y##*/}.raw.bcf $RESULTSMPILE/${y##*/}.raw.bcf
#			else
#				echo "${y}/${y##*/}.raw.bcf does not exist."
#			fi
#			if [ -f ${y}/${y##*/}_merged.muts.vcf ];then
#				if [ ! -f ${2}/${y##*/}_merged.muts.vcf ];then
#					cp ${y}/${y##*/}_merged.muts.vcf ${2}/
#				else
#					echo "File has already been copied"
#				fi
#			else
#				echo "${y##*/}_merged.muts.vcf does not exist."
#			fi

			if [ -f  $y/run1/${y##*/}_vs_DDv3_1is.snps.ids.vcf.gz ];then
			  cp $y/run1/${y##*/}_vs_DDv3_1is.snps.ids.vcf.gz ${RESULTSSNPS} 
			else

				echo "${y##*/}__vs_DDv3_1is.snps.ids.vcf.gz does not exist"

			fi
                done
}
function move_pindel_results()
{
for z in $(ls -d "$PINDELCGP"/*-Dog)
	do 
		if [ -f ${PINDELCGP}/${z##*/}/${z##*/}_vs_DDv3_1is.flagged.vcf ];then
			if [ ! -f ${RESULTSPINDEL}/${z##*/}_vs_DDv3_1is.flagged.vcf ];then
				echo "Copying ${z##*/}_vs_DDv3_1is.flagged.vcf to $RESULTSPINDEL"
				cp "$PINDELCGP"/${z##*/}/${z##*/}_vs_DDv3_1is.flagged.vcf "$RESULTSPINDEL"
			else
				echo "File has already been copied."
			fi
		else
			echo "${z##*/}_vs_DDv3_1is.flagged.vcf does not exist."
		fi
	done
}
function polyc_expansion()
{
	cd ${POLYC}/
	for q in $(ls -d "$TUMOUR"/* "$HOST"/*);do
		sample_name=$(echo ${q##*/} | cut -d'.' -f1)
		#check that bam file exists. important.
		if [ -f ${q}/$sample_name.sorted.rmdup.bam ];then
			#extract mapped reads with an unmapped mate
			if [ ! -f ${POLYC}/${sample_name}.dat ];then
				echo "Making data file for $sample_name"
				samtools view -b -f8 -F4 $q/$sample_name.sorted.rmdup.bam > $q/$sample_name.polyc.bam
				samtools depth -b $runHOME/bed/MT.bed ${q}/$sample_name.polyc.bam | awk 'BEGIN { prev_chr="";prev_pos=0;} { if($1==prev_chr && prev_pos+1!=int($2)) {for(i=prev_pos+1;i<int($2);++i) {printf("%s\t%d\t0\n",$1,i)
;}} print; prev_chr=$1;prev_pos=int($2);}' | cut -f2,3 > ${POLYC}/$sample_name.dat
			else
                                echo "$sample_name.dat already exists"
                        fi
			# if CIGAR indicates that read is soft clipped or hard clipped then extract. extract header lines also.	
			if [ ! -f ${POLYC}/${sample_name}_clipped.dat ];then
				echo "Extracting clipped reads from $sample_name.sorted.rmdup.bam"
				samtools view -h $q/$sample_name.sorted.rmdup.bam | awk '{if ($1 ~ /@/ || $6 ~ /S/ || $6 ~ /H/) print $0 }' > $q/${sample_name}_clipped.sam
				echo "Done!"
				echo "Converting clipped reads sam to bam."
				samtools view -bS $q/${sample_name}_clipped.sam > $q/${sample_name}_clipped.bam
				echo "Done!"
				#awk part fills in the gaps where 0 reads aligned to that base
				samtools depth -b $runHOME/bed/MT.bed ${q}/${sample_name}_clipped.bam | awk 'BEGIN { prev_chr="";prev_pos=0;} { if($1==prev_chr && prev_pos+1!=int($2)) {for(i=prev_pos+1;i<int($2);++i) {printf("%s\t%d\t0\n",$1,i);}} print; prev_chr=$1;prev_pos=int($2);}' | cut -f2,3 > ${POLYC}/${sample_name}_clipped.dat
			else
				echo "${sample_name}_clipped.dat already exists."
			fi
			#plot position vs read depth
			if [ ! -f ${POLYC}/$sample_name.png ];then
				echo "Making plot for $sample_name"
				gnuplot -e "samplename='$sample_name';filename='$sample_name.dat'" polyc.gnu > ${POLYC}/$sample_name.png 
			else 
				echo "$sample_name.png already exists"
			fi
			if [ ! -f ${POLYC}/${sample_name}_clipped.png ];then
                                echo "Making plot for $sample_name"
                                gnuplot -e "samplename='$sample_name';filename='${sample_name}_clipped.dat'" polyc.gnu > ${POLYC}/${sample_name}_clipped.png  
                        else 
                                echo "${sample_name}_clipped.png already exists"
                        fi
		else
			echo "$sample_name.sorted.rmdup.bam does not exist"
			
		fi	
	done

	# convert ${POLYC}/*T*png $runHOME/polyc_graphs/${today}_tumour_polyc_graphs.pdf
	# convert ${POLYC}/*H*png $runHOME/polyc_graphs/${today}_host_polyc_graphs.pdf

}

function bam_coverage_stats(){
	echo "Computing average coverage for all bam files."
	for r in $(ls -d "$TUMOUR"/* "$HOST"/*);do
		samplename=$(echo ${r##*/} | cut -d'.' -f1)
		echo $samplename >> ${today}_average_coverage.txt
		samtools depth -b ${BED}/MT.bed ${r}/$samplename.bam | awk '{sum+=$3; $4=sum/16726}{print $4}' | tail -n 1 >> ${runHOME}/${today}_average_coverage.txt
	done
	echo "Done!"
}
function copy_to_dogndevil(){
	for s in $(ls -d "$TUMOUR"/* "$HOST"/*);do
                samplename=$(echo ${s##*/} | cut -d'.' -f1)
		ls -l $s/${samplename}.bam
		cp $s/${samplename}.bam /nfs/dog_n_devil/dog/mitochondria/bam
		cp $s/${samplename}.bam.bai /nfs/dog_n_devil/dog/mitochondria/bam
		ls -l /nfs/dog_n_devil/dog/mitochondria/bam/${samplename}.bam
	done

}

canapps_download T $TUMOUR
canapps_download H $HOST
check_file_size
index_sort $TUMOUR/*
check_coverage 
list_unmatched
pindel_bams
run_pindelcgp
make_pindel_bedfile
run_cavemanwrapper ${TUMOUR} $TUMOUR no_analysis_tumour.bed
run_cavemanwrapper ${HOST} $HOST no_analysis_host.bed
run_mpileup
move_caveman_results "$HOST" "$RESULTSMUTS"/all
move_caveman_results "$TUMOUR" "$RESULTSMUTS"/all
move_pindel_results
polyc_expansion
bam_coverage_stats
copy_to_dogndevil
