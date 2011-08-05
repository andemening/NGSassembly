#!/bin/bash

# Functions file
# build basespace binary from genome refseq
# build colorspace binary from genome refseq
# create jump database from colorspace binary
# align genome in colorspace 
# sorting of reads according to coordinates, resolve mate-pairs
# convert alignment file to BAM 
# mark duplicates with Picard MarkDuplicates
# merge all alignments and pipe to FreeBayes to call SNPs

################# Parallel threads ####################
function Build {					## Create binary file in basespace from genome FASTA file
		# -fr,-fr2 bases for 1st mate and 2nd mate
		# -fq,-fq2 base qualities for 1st mate and 2nd mate
		# -assignQual assign a quality value to all bases
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
		
		# after quality filtering >20 and mate-pair pairing, 163.00 reads passed out of 338,21 reads total
		# ~/mosaik-aligner/bin/MosaikBuild -fr gallus_HS_MP_2x50bp_T_mp_F3.csfasta -fq gallus_HS_MP_2x50bp_T_mp_F3_QV.qual -fr2 gallus_HS_MP_2x50bp_T_mp_R3.csfasta -fq2 gallus_HS_MP_2x50bp_T_mp_R3_QV.qual -st solid -ln HIGH_LINE -id HS_MP_2x50bp -ds "gallus high HS 2x50bp mate-pairs" -sn "gallus gallus domesticus" -mfl 400 -out gallus_HS_2x50bp.dat | tee -a build.gallus_HS_2x50bp.log

	echo -e "\tCreating binary file in basespace for $project_name..."		
	$mosaik_dir/MosaikBuild -fr $refseq_dir/$file -oa $results_dir/$file".dat" | tee -a $results_dir/build.log
	$mosaik_dir/fixlog $results_dir/build.log
}


function BuildCS {    				## Create binary file in colorspace from genome FASTA file
	echo -e "\tCreating binary file in colorspace for $project_name..."
	$mosaik_dir/MosaikBuild -fr $refseq_dir/$file -oa $results_dir/$file".cs.dat" -cs | tee -a $results_dir/buildCS.log
	$mosaik_dir/fixlog $results_dir/buildCS.log
}


function JumpDB {					## Create jump database for genome 
		# -kd, keeps keys database on disk
		# -mem, amount of RAM used when sorting, hard-coded to 70 GB
		# -mhp 100, maximum hash positions, default 100
		
	st=$SECONDS
	echo -e "\tCreating jumps database for $project_name..."
	$mosaik_dir/MosaikJump -ia $results_dir/$file".cs.dat" -hs $hash_size -mhp $mhp -out $results_dir/$file".jmp" | tee -a $results_dir/mosaik.log
	$mosaik_dir/fixlog $results_dir/mosaik.log
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tJumps database done in $exectime seconds." | tee -a $main_dir/pipeline.log
}


function Align {					## Align
		# -in SOLiD reads dat file; -ia refseq in colorspace; -ibs -ibs refseq in basespace
		# -ls search radius for mate-pairs
		# -bw banded Smith-Waterman alignment algorithm 
		# -mhp maximum hash positions 
		# -act alignment candidate threshold
		# -mm mismatches
		# -mmp threshold of ALIGNMENT QUALITY, AlignmentQualityThreshold, number between 0 & 99, =10 for conservative filtering
		# -j jumpdatabase, $results_dir/$chromosome/$chromosome".jmp
 		# 
	

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
		
	st=$SECONDS
	echo -e "\tAligning genome $genome..."
	
	read_length=`grep -q "35bp" $genome > /dev/null; echo $?`
	if [ "$read_length" -eq "0" ]; then
		mismatches=$SE_mismatches
		bandwidth=$SE_bandwidth
		act=$SE_act
	else
		mismatches=$MP_mismatches
		bandwidth=$MP_bandwidth
		act=$MP_act
	fi

	$mosaik_dir/MosaikAligner -in $reads_dir/$genome -ia $results_dir/$file".cs.dat" -ibs $results_dir/$file".dat" -hs $hash_size -mm $mismatches -bw $bandwidth -act $act -mhp $mhp -m all -a all -p $proc -out $results_dir/$genome/$file".aligned.dat" -j $results_dir/$file".jmp" | tee -a $results_dir/mosaik.log
	$mosaik_dir/fixlog $results_dir/mosaik.log
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tAlignment done in $exectime seconds." | tee -a $main_dir/pipeline.log
}


function Sort {					## Sorting, two different approaches for SE and MP
		# simply sorts SE
		# resolves MP, fragment length is calculated using confidence interval, min,max,mean of distribution is returned
		# -mem number of alingments in memory, sets sorting cache size, default 6000000
		# -nu include non-unique reads, ONLY for Single-End
		
		# -afl allow all fragment lengths for PE/MP reads
		# -sa sample all unique vs non-unique fragment lengths
		# -rmm resolve multiple vs multiple read pairs
		# -dup <Directory>, filters out duplicate reads, DISCARDS reads with lower alignment quality score
		
	st=$SECONDS
	echo -e "\tSorting aligned reads for $genome..."
	$mosaik_dir/MosaikSort -in $results_dir/$genome/$file".aligned.dat" -mem 8000000 -nu -afl -sa -rmm -out $results_dir/$genome/$file".aligned.sorted.dat" | tee -a $results_dir/mosaik.log
	$mosaik_dir/fixlog $results_dir/mosaik.log
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tSorting done in $exectime seconds." | tee -a $main_dir/pipeline.log
}


function Text {					## Convert assembly archive to BAM file format for analysis in FreeBayes
	st=$SECONDS
	eco -e "\tConverting alignment file for $genome to BAM file…"
	$mosaik_dir/MosaikText -in $results_dir$genome/$file".aligned.sorted.dat" -bam $results_dir/$genome/$file".aligned.sorted.bam" | tee -a $results_dir/mosaik.log
	$mosaik_dir/fixlog $results_dir/mosaik.log
	en=$SECONDS
	exectime=$((en -st))
	echo -e "\tConversion done in $exectime seconds." | tee -a $main_dir/pipeline.log
}


# Index duplicate reads in database contained in folder (flag -od)
function MarkDuplicates {			## Remove duplicate PCR products to improve SNP detection
		# -od, output directory for duplicate database
		
	st=$SECONDS
	echo -e "\tMarking duplicate reads for $genome..."
		# $mosaik_dir/MosaikDupSnoop -in $results_dir/$chromosome/$genome/$chromosome".aligned.dat" -afl -rmm -od $results_dir/$chromosome/"Duplicates"
		
	java -d64 -Xmx18G -jar $mosaik_dir/MarkDuplicates.jag INPUT=$results_dir/$genome/$file".aligned.sorted.bam" OUTPUT=$results_dir/$genome/$file".aligned.sorted.Marked.bam" METRICS=$results_dir/$genome/$file".Marked.metrics.txt" REMOVE_DUPLICATES=false ASSUME_SORTED=true | tee -a $results_dir/$genome/markduplicates.log
	en=$SECONDS
	exectime=$((en - st))
	echo -e "\tMarked duplicates in $exectime seconds." | tee -a $main_dir/
}


#function Merge {					## Merge all sequencing lines (high line and low line) for overall assembly coverage 
#		# -mem number of alignments in memory, sets sorting cache, default 6000000
#
#	st=$SECONDS
#	echo -e "\t Merging all lines for $chromosome…"
#
#	dedup.sorted_files=`find $results/$chromosome -name \*.dedup.sorted.dat -print`	# lists all bam-files in current chromosome
#	infiles=""
#	for merging_name in $dedup.sorted_files; do
#		infiles="$infiles-in $merging_name "							# generate in variable for merging command
#	done
#
#	$mosaik_dir/MosaikMerge $infiles -mem 12000000 -out $results_dir/$chromosome/$chromosome".merged.dat" | tee -a $results_dir/merge.coverage.log
#	$mosaik_dir/fixlog $results_dir/merge.coverage.log
#	en=$SECONDS
#	exectime=$((en - st))
#	echo -e "\tMerge done in $exectime seconds." | tee -a $main_dir/pipeline.log
#}


#function Coverage {				## Calculate coverage for all sequencing lines (high line and low line)
#		# -in alignment file
#		# -ia reference sequence colorspace dat file
#		# -u limits coverage analysis to unique reads - duplicate reads already removed
#		# -mc minimum coverage, default 1
#		# -cg plots graph using gnuplot, also containing pdf - takes time for large chromosomes
#		# -od output directory for analysis
#		# -ngc skips making graphs and coverage files
#		
#	st=$SECONDS
#	echo -e "\tCalculating coverage for all lines of chromosome $chromosome (unique reads only)..."
#	$mosaik_dir/MosaikCoverage -in $results_dir/$chromosome/$chromosome".merged.dat" -ia $refseq_dir/$chromosome/$chromosome".dat" -od $results_dir/$chromosome/"coverage" | tee -a $results_dir/merge.coverage.log
##	~/mosaik-aligner/bin/MosaikCoverage -ia ~/chr2/chr2.cs.dat -in chr2.aligned.sorted.dat -u -mc 3 -cg -od coverage
#	$mosaik_dir/fixlog $results_dir/merge.coverage.log
#	en=$SECONDS
#	exectime=$((en - st))
#	echo -e "\tPlotting coverage done in $exectime seconds." | tee -a $main_dir/pipeline.log
#}


#function Assemble {				## Multiple sequence alignment and assembly 
#		# -in alignment file
#		# -ia reference sequence colorspace dat file
#		# -out filename stub/prefix for output file
#		# -f output file format, default ACE
#	
#	st=$SECONDS
#	echo -e "\tAssembling chromosome $chromosome in line $genome..."
#	$mosaik_dir/MosaikAssembler -in $results_dir/$chromosome/$genome/$chromosome".sorted.dedup.dat" -ia $results_dir/$chromosome/$chromosome".cs.dat" -f ace -out assembly_ | tee -a $results_dir/$chromosome/mosaik.log
#	$mosaik_dir/fixlog $results_dir/$chromosome/mosaik.log
#	en=$SECONDS
#	exectime=$((en - st))
#	echo -e "\tAssembly done in $exectime seconds." | tee -a $main_dir/pipeline.log
#}


function STM { 				## These functions run in sequential fashion
#			Merge					# Merge all chromosomes in each line to get complete coverage; total coverage for all pooled sequences (important for SNP statistics?) 
#			Coverage				# Calculate coverage
#		Assemble					# Assemble alignment chromosome by chromosome and output ACE file
		Sort						# Sort SE reads; sort and resolve MP reads
		Text						# Convert binary alignment DAT to BAM
		MarkDuplicates				# Mark PCR duplicates
}


####################### SNP call pipe #####################
# all previous functions must have completed before this step

function CallSNP {				## Calling SNPs and short INDELS
	SNPstart=$SECONDS
	echo -e "\tSNP calling in $genome..."

#	# Add read group names with bamaddrg
#	# This is not required in current analysis of 3 high-line, 3 low-line
#
#	# bamaddrg settings: 
#		# -s sample name, set to chromosome
#		# -r read group, set to genome/reads filename (gallus_high_35_1, gallus_low_50, ...)
#	bam_files=`find $results/$chromosome -name \*.bam -print`	# lists all bam-files in current chromosome
#	infiles=""
#	for sample_name in $bam_files; do
#		$mosaik_dir/bamaddrg -b $sample_name -s $chromosome -r $genome				# appends each $genome as read group name for for corresponding bam-files
#		# infiles is string with all input files for bayesian analysis
#		# $infiles = "-in /complete/path/file1.bam -in /complete/path/file2.bam -in ... "
#		infiles="$infiles-in $sample_name "									# generate list of bam-files for 
#	done

	
	bam_files=`find $results -name \*.Marked.bam -print | awk -F "/" '{print $NF}'`	# lists all bam-files in current line
	infiles=""
	for merging_name in $bam_files; do
		$mosaik_dir/bamtools index $results_dir/$merging_name
		infiles="$infiles-in $merging_name "								# generate in variable for merging command
	done

	# Pipe entire SNP call command from stdout to stdin; bamtools merge -> bamleftalign -> freebayes
	# Bamtools merge piped to bamleftalign for homogenization of positional distribution of IN/DELs, computationally inexpensive
	# bamtools merge xxx | bamleftalign xxx | freebayes --stdin -f $refseq_dir/$chomosome/$file --ploidy diploid --indel -dd --pvar $PSL --min-coverage $CAL -v $results_dir/$chromosome/$chromosome".called.vcf" --log $results_dir/$chromosome/$genome$/freebayes.merged.leftaligned.log

	# bamleftalign only required for MP reads when INDEL calling
#	$mosaik_dir/bamleftalign -f $refseq_dir/$file/


	# no sample names are given -> using name of BAM file as sample name
	$mosaik_dir/bamtools merge $infiles | $mosaik_dir/freebayes --stdin -f $refseq_dir/$file --pooled --ploidy 22 -dd --pvar $PSL --min-coverage $CAL --min-mapping-quality $MAQ --ignore-reference-allele --use-mapping-quality -v $results_dir/$line/$file".called.vcf" --log $results_dir/$line".called.log"
		# -f fasta refseq, $file is current chromosome basename in iteration 		
		# -v, output file
		# -log, log filename
		# --ploidy, default diploid
		# -d, debugging output
		# --pvar, reporting probaility limit
		# --min-coverage, least required coverage to process site
		# --min-mapping-quality, minimum default 30
		# --use-mapping-quality
		# --ignore-reference-allele 

		
	SNPend=$SECONDS
	SNPtime=$((SNPend - SNPstart))
	echo -e "\tSNP calling done in $exectime seconds." | tee -a $main_dir/pipeline.log
}


function CallINDEL {				## Calling SNPs and short INDELS
	SNPstart=$SECONDS
	echo -e "\tINDEL calling in $genome..."

	bam_files=`find $results/$genome -name \*.Marked.bam -print`	# lists all bam-files in current line
	infiles=""
	for merging_name in $bam_files; do
		infiles="$infiles-in $merging_name "							# generate in variable for merging command
	done

	# bamleftalign only required for MP reads when INDEL calling
	$mosaik_dir/bamleftalign -f $refseq_dir/$genome/

	$mosaik_dir/bamtools index $results_dir/$file
	$mosaik_dir/bamtools merge $infiles | $mosaik_dir/freebayes --stdin -f $refseq_dir/$file".fa" --ploidy 2 -dd --pvar $PSL --min-coverage $CAL --min-mapping-quality $MAQ --ignore-reference-allele --use-mapping-quality -v $results_dir/$genome/$file".called.vcf" --log $results_dir/$genome
		# -f fasta refseq, $file is current chromosome basename in iteration 		
		# -v, output file
		# -log, log filename
		# --ploidy, default diploid
		# -d, debugging output
		# --pvar, reporting probaility limit
		# --min-coverage, least required coverage to process site
		# --min-mapping-quality, minimum default 30
		# --use-mapping-quality
		# --ignore-reference-allele 

		
	SNPend=$SECONDS
	SNPtime=$((SNPend - SNPstart))
	echo -e "\tINDEL calling done in $exectime seconds." | tee -a $main_dir/pipeline.log
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
			regeneratequeue					# new queue if some job/PID has finished
			break							# exit if-loop
		fi
	done
}

