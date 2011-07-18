#!/bin/bash

## Fix log paths, each analysis log in separate matching folder
## Fix input output filenames
## Do RemoveDuplicates, Merge and Text functions

################# Parallel threads ####################
function Build {					## Create binary file in basespace

	echo -e "\tCreating binary file in basespace..."
		# -fr,-fr2 bases for 1st mate and 2nd mate
		# -fq,-fq2 base qualities for 1st mate and 2nd mate
		# -assignQual assign
		# -oa outfile for reference sequence; -out outfile when converting reads
		# -st specifies sequencing technology; -ds description; -ln library name
		# -mfl specify median fragment length for mate-pairs

		# -ga genome assembly ID, eg. galGal3
		# -sn species name, eg. "gallus gallus domesticus"
		# -sam sample name, HS_MP_2x50bp 	MW_MP_2x50bp 	HS_SE_35bp_1 	HS_SE35bp_2 	MW_SE_35bp_1 	MW_SE_35bp_2
		
		# Archive gallus_high_2x50bp_matepairs_csfasta.qual.tar.gz contains: 
		# ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_F3.csfasta
		# ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_F3_QV.qual
		# ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_R3.csfasta
		# ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_R3_QV.qual
		# Library: HIGH_LINE	LOW_LINE
		# Group: HS_MP_2x50bp 	MW_MP_2x50bp 	HS_SE_35bp_1 	HS_SE35bp_2 	MW_SE_35bp_1 	MW_SE_35bp_2

		# -assignQual <value> = assign <value> quality to every base
		# ~/mosaik-aligner/bin/MosaikBuild -fr ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_F3.csfasta.gz -fq ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_F3_QV.qual.gz -fr2 ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_R3.csfasta.gz -fq2 ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_R3_QV.qual.gz -st solid -ln HIGH_LINE -id HS_MP_2x50bp -ds "gallus high HS 2x50bp mate-pairs" -sn "gallus gallus domesticus" -mfl 400 -out gallus_HS_2x50bp.dat | tee -a build.gallus_HS_2x50bp.log
		#  -ln HIGH_LINE -id HS_MP_2x50bp -ds "gallus high HS 2x50bp mate-pairs" -sn "gallus gallus domesticus"

		# using old version MosaikBuild 1.0.1388
		# ~kiero/gallus/bin/MosaikBuild -fr ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_F3.csfasta -fq ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_F3_QV.qual -fr2 ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_R3.csfasta -fq2 ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_R3_QV.qual -st solid -ln HIGH_LINE -id HS_MP_2x50bp -ds "gallus high HS 2x50bp mate-pairs" -sn "gallus gallus domesticus" -mfl 400 -out Mosaik1.0_gallus_HS_2x50bp.dat | tee -a build.Mosaik1.0_gallus_HS_2x50bp.log

		# Quality filter and mate-pair pairing
		# mean Phred quality score 20 per read
		# perl ~/SOLiD_preprocess_meanFilter_v2.pl -i mp -f ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_F3.csfasta -g ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_F3_QV.qual -r ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_R3.csfasta -s ugc_208_209_MP_2x50bp_solid0105_20091127_ugc_208_R3_QV.qual -p on -q 20 -a on -v on -o gallus_HS_MP_2x50bp | tee -a preprocess_meanFilter.log
		
		# after quality filtering >20 and mate-pair pairing, 163.00 Mb passed out of 338,21 Mb
		# ~/mosaik-aligner/bin/MosaikBuild -fr gallus_HS_MP_2x50bp_T_mp_F3.csfasta -fq gallus_HS_MP_2x50bp_T_mp_F3_QV.qual -fr2 gallus_HS_MP_2x50bp_T_mp_R3.csfasta -fq2 gallus_HS_MP_2x50bp_T_mp_R3_QV.qual -st solid -ln HIGH_LINE -id HS_MP_2x50bp -ds "gallus high HS 2x50bp mate-pairs" -sn "gallus gallus domesticus" -mfl 400 -out gallus_HS_2x50bp.dat | tee -a build.gallus_HS_2x50bp.log

		
	$mosaik_dir/MosaikBuild -fr $refseq_dir/$file -oa $results_dir/$chromosome/$chromosome".dat" -mfl $mfl | tee -a $results_dir/$chromosome/
	$mosaik_dir/fixlog $results_dir/$chromosome/
}


function BuildCS {    				## Create binary file in colorspace
	echo -e "\tCreating binary file in colorspace..."
	$mosaik_dir/MosaikBuild -fr $refseq_dir/$file -oa $results_dir/$chromosome/$chromosome".cs.dat" -mfl $mfl -cs | tee -a $results_dir/$chromosome/
	$mosaik_dir/fixlog $results_dir/$chromosome/
}


function JumpDB {					## Create jump database
	st=$SECONDS
	echo -e "\tCreating jumps database..."
		# -kd, keeps keys database on disk
		# -mem, amount of RAM used when sorting, hard-coded to 70 GB

	$mosaik_dir/MosaikJump -ia $results_dir/$chromosome/$chromosome".cs.dat" -hs $hash_size -mem 70 -out $results_dir/$chromosome/$chromosome".jmp" | tee -a $results_dir/$chromosome/
	$mosaik_dir/fixlog $results_dir/$chromosome/
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tJumps database done in $exectime seconds." | tee -a $main_dir/
}


function Align {					## Align
	st=$SECONDS
	echo -e "\tAligning genome $genome..."
		# -in SOLiD reads in; -ia refseq in colorspace; -ibs -ibs refseq in basespace

		# -ls search radius for mate-pairs
		# -bw banded Smith-Waterman alignment algorithm 
		# -mhp maximum hash positions 
		# -act alignment candidate threshold
		# -mm mismatches
		#-j $results_dir/$chromosome/$chromosome".jmp
	
		# single end, practical example
		# ~/mosaik-aligner/bin/MosaikAligner -in /home/gallus/gallus_high_35_1.dat -ia chr2.cs.dat -ibs chr2.dat -a all -m all -hs 15 -p 16 -bw 13 -act 20 -mm 4 -mhp 100 -j chr2.cs.jmp -out results.gal_high_35bp_1.16cores/chr2.aligned.dat | tee -a results.gal_high_35bp_1.16cores/aligner.log
		
		# mate-pairs, practical example
		# ~/mosaik-aligner/bin/MosaikAligner -in /home/gallus/gallus_high_50.dat -ia chr2.cs.dat -ibs chr2.dat -a all -m all -hs 15 -p 16 -bw 17 -act 25 -mm 6 -ls 100 -mhp 100 -j chr2.cs.jmp -out results.gal_high_50bp.16cores/chr2.aligned.dat | tee -a results.gal_high_50bp.16cores/aligner.log

		# Alignment after quality filtering >20 and mate-pair pairing (163.00 Mb passed out of 338,21 Mb)
		# ERROR: This sequencing technology is not currently supported for local alignment search.
		# ~/mosaik-aligner/bin/MosaikAligner -in /home/gallus/gallus_HS_2x50bp.dat -ia chr2.cs.dat -ibs chr2.dat -a all -m all -hs 15 -p 16 -bw 17 -act 25 -mm 6 -ls 100 -mhp 100 -j chr2.cs.jmp -out results.gallus_HS_2x50bp/chr2.aligned.dat | tee -a results.gallus_HS_2x50bp/aligner.log
	
		# Alignment of Mosaik1.0_gallus_HS_2x50bp built using MosaikBuild 1.0.1388 from F3.csfasta + F3.qual and R3.csfasta + R3.qual 
		# ERROR: This sequencing technology is not currently supported for local alignment search.
		# Working after omitting flag -ls 100 (radius search of median fragment length), local alignment search of SOLiD mate-pairs in colorspace not supported ?? 
		# ~/mosaik-aligner/bin/MosaikAligner -in /home/gallus/Mosaik1.0_gallus_HS_2x50bp.dat -ia chr2.cs.dat -ibs chr2.dat -a all -m all -hs 15 -p 16 -act 25 -mm 6 -bw 17 -ls 100 -mhp 100 -j chr2.cs.jmp -out results.Mosaik1.0_gallus_HS_2x50bp/chr2.aligned.dat | tee -a results.Mosaik1.0_gallus_HS_2x50bp/aligner.log
		
		# Alignment of Mosaik1.0_gallus HS_2x50bp.mfl3528 built using MosaikBuild 1.0.1388 with flag -mfl 3528 (mean median fragment length) 
		# ~/mosaik-aligner/bin/MosaikAligner -in /home/gallus/Mosaik1.0_gallus_HS_2x50bp.mfl3528.dat -ia chr2.cs.dat -ibs chr2.dat -a all -m all -hs 15 -p 16 -act 25 -mm 6 -bw 17 -ls 600 -mhp 100 -j chr2.cs.jmp -out results.Mosaik1.0_gallus_HS_2x50bp.mfl3528/chr2.aligned.dat | tee -a results.Mosaik1.0_gallus_HS_2x50bp.mfl3528/aligner.log
		
	$mosaik_dir/MosaikAligner -in $reads_dir/$genome -ia $results_dir/$chromosome/$chromosome".cs.dat" -ibs $results_dir/$chromosome/$chromosome".dat" -hs $hash_size -mm $max_mismatches -act $algn_thresh -m all -a all -p $proc -ls $radius -out $results_dir/$chromosome/$genome/$chromosome".aligned.dat" | tee -a $results_dir/$chromosome/
	$mosaik_dir/fixlog $results_dir/$chromosome/
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tAlignment done in $exectime seconds." | tee -a $main_dir/
}


function Sort {					## Sorting, two different approaches for SE and MP
	# simply sorts SE
	# resolves MP, fragment length is calculated using confidence interval, min,max,mean of distribution is returned
	st=$SECONDS
	echo -e "\tSorting aligned reads for $genome..."
	$mosaik_dir/MosaikSort -in $results_dir/$chromosome/$genome/$chromosome".aligned.dat" -out $results_dir/$chromosome/$genome/$chromosome".sorted.dat" | tee -a $results_dir/$chromosome/
	
	# 16 processes at same time
	# ~/mosaik-aligner/bin/MosaikSort -in chr2.aligned.dat -mem 25000000 -afl -sa -out chr2.aligned.sorted.dat | tee -a mosaik.log
	
	$mosaik_dir/fixlog $results_dir/$chromosome/mosaik.log
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tSorting done in $exectime seconds." | tee -a $main_dir/
}



function RemoveDuplicates { 	## Remove duplicate PCR products to improve SNP detection
	$mosaik_dir/MosaikDupSnoop -in $results_dir/$chromosome/$genome/$chromosome".sorted.dat" -out $results_dir/$chromosome/$chromosome".sorted.dedup.dat" | tee -a $results_dir/$chromosome/
}


function Merge {					## Merge all sequencing lines (high line and low line) for overall assembly coverage 
	st=$SECONDS
	echo -e "\t Merging all lines for $chromosome…"

	$mosaik_dir/MosaikMerge 
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tMerge done in $exectime seconds." | tee -a $main_dir/
}


function Coverage {				## Calculate coverage
	st=$SECONDS
	echo -e "\tCalculating coverage for whole genome (unique reads only)..."
		# -u limits coverage analysis to unique reads - duplicate reads already removed
		# -mc minimum coverage, default 1
		# -cg plots graph using gnuplot, also containing pdf - takes time for large chromosomes
		# -od output directory for analysis
		# -ngc skips makign graphs and coverage files
		
	$mosaik_dir/MosaikCoverage -in $results_dir/$chromosome/$genome/$chromosome".sorted.dedup.dat" -ia $results_dir/$chromosome/$chromosome".dat" -od $results_dir/$chromosome/$genome/"coverage" | tee -a $results_dir/$chromosome/

		# practical example
		# ~/mosaik-aligner/bin/MosaikCoverage -ia ~/chr2/chr2.cs.dat -in chr2.aligned.sorted.dat -u -mc 3 -cg -od coverage

	$mosaik_dir/fixlog $results_dir/$chromosome/
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tPlotting coverage done in $exectime seconds." | tee -a $main_dir/
}


function Assemble {				## Assembling
	st=$SECONDS
	echo -e "\tAssembling $genome..."
		# -cf creates FASTA read file for GigaBayes, old version
	
	$mosaik_dir/MosaikAssembler -in $results_dir/$chromosome/$genome/$chromosome".sorted.dedup.dat" -ia $results_dir/$chromosome/$chromosome".cs.dat" -f gig -out $results_dir/$chromosome/$genome/"assembly_"$chromosome".ace" | tee -a $results_dir/$chromosome/
	$mosaik_dir/fixlog $results_dir/$chromosome/
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tAssembly done in $exectime seconds." | tee -a $main_dir/
}


function Text {					## Convert assembly archive to BAM file format for analysis in FreeBayes
	$mosaik_dir/MosaikText -in $results_dir/$chromosome/$genome/"assembly_"$chromosome".ace" -bam $results_dir/$chromosome/$genome/$chromosome".assembly.bam" | tee -a 
}

################# All threads at once ####################
# all previous functions must have completed before this step

function CallSNP {				## Calling SNPs and short INDELS
	SNPstart=$SECONDS
	echo -e "\tSNP/INDEL calling in $genome..."

	bam_files=`find $results/$chromosome -name \*.bam -print`	# lists all bam-files in current chromosome
	

	# Add read group names with bamaddrg
	# bamaddrg settings: 
		# -s sample name, set to chromosome
		# -r read group, set to genome/reads filename (gallus_high_35_1, gallus_low_50, ...)
	infiles=""
	for sample_name in $bam_files; do
		$mosaik_dir/bamaddrg -b $sample_name -s $chromosome -r $genome			# appends each $genome as read group name for for corresponding bam-files
		# infiles is string with all input files for bayesian analysis
		# $infiles = "-in /complete/path/file1.bam -in /complete/path/file2.bam -in ... "
		infiles="$infiles-in $sample_name "
	done
	
	# Bamtools merge piped to bamleftalign for homogenization of positional distribution of IN/DELs, computationally inexpensive
	# Pipe entire SNP call command from stdout to stdin; bamtools merge -> bamleftalign -> freebayes
	# bamtools merge xxx | bamleftalign xxx | freebayes --stdin -f $refseq_dir/$chomosome/$file --ploidy diploid --indel -dd --pvar $PSL --min-coverage $CAL -v $results_dir/$chromosome/$chromosome".called.vcf" --log $results_dir/$chromosome/$genome$/freebayes.merged.leftaligned.log

	$mosaik_dir/bamtools merge$infiles | $mosaik_dir/bamleftalign -f $refseq_dir/$chomosome/$file | $mosaik_dir/freebayes --stdin -f $refseq_dir/$chomosome/$file 

	$mosaik_dir/freebayes -b $results_dir/$chromosome/$genome/$bam_file -f $refseq_dir/$chomosome/$file --ploidy diploid --indel -mnps --dd --pvar $PSL --min-coverage $CAL -v $results_dir/$chromosome/$chromosome".called.vcf" --log $results_dir/$chromosome/$genome/freebayes.log
	SNPend=$SECONDS
	SNPtime=$((SNPend - SNPstart))
	echo -e "\tSNP/INDEL calling done in $exectime seconds." | tee -a $main_dir/
}




################# Multi-core functions ####################
function queue {
	QUEUE="$QUEUE $1"
	NUM=$(($NUM+1))
}


function regeneratequeue {
	OLDREQUEUE=$QUEUE
	QUEUE=""
	NUM=0
	for PID in $OLDREQUEUE
	do
		if [ -d /proc/$PID ] ; then				# If directory exists = PID running
			QUEUE="$QUEUE $PID"
			NUM=$(($NUM+1))
		fi
	done
}


function checkqueue {
	OLDCHQUEUE=$QUEUE
	for PID in $OLDCHQUEUE
	do
		if [ ! -d /proc/$PID ] ; then 			# If directory does not exist -> start new process
			regeneratequeue
			break
		fi
	done
}

