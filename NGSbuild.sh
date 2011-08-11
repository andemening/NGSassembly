#!/bin/bash

####################
# Build binary Mosaik dat files containing filtered raw reads from SOLiD sequencing

# Rutger's perl script for quality filtering, keep =< 10 Phred+33
# analyze QV distribution, observe the difference between base quality and mapping quality

# Read groups, sample names need to be set in MosaikBuild

	bam_files=`find $results -name \*.Marked.bam -print | awk -F "/" '{print $NF}'`	# lists all bam-files in current line
	infiles=""
	for merging_name in $bam_files; do
		case $merging_name in 
			RG_basename=`echo $merging_name | awk -F "/" '{print $NF}' | awk -F "." '{print $1}'`
			*"high"*)
				infiles="$infiles -b $merging_name -s high -r $RG_basename"
			*"low"*)
				infiles="$infiles -b $merging_name -s $low -r $RG_basename"
		esac
		$mosaik_dir/bamtools index $results_dir/$merging_name
	done
	$mosaik_dir/bamaddrg $infiles | $mosaik_dir/freebayes --stdin -f $refseq_dir/$file --pooled --ploidy $ploidy -dd --pvar $poly_prob --min-coverage $CAL --min-mapping-quality $MAQ --ignore-reference-allele --use-mapping-quality -v $results_dir/$line/$line"_"$file".called.vcf" --log $results_dir/$line".called.log"
