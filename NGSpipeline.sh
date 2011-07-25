#!/bin/bash

###########################################################
# Genome assembly pipeline using Mosaik suite of programs #
# and calling SNPs and short INDELs using FreeBayes       #
# Automated and interactive bash script using multiple    #
# threads, analysis done chromosome by chromosome         #
# Script by Marcin Kierczak                               #
# Marcin.Kierczak@hgen.slu.se, date: 18.05.2010           #
# Modified and expanded by: Andreas.E.Lundberg@gmail.com  #
# Date: 14.07.2011                                        #
###########################################################

# Genome assembly pipeline using Mosaik suite of programs and SNP calling with FreeBayes
# Automated and interactive bash script using multiple threads for speed
#
# Program assembles one chromosome at a time from individual chromosome fasta files. \n
# Generates individual chromosome directory in results directory.
# First set up working directories and pipeline variables !!
# Alignment is done in linear processing fashion using multiple cores.
# All other steps are run IN PARALLEL utilizing user selectible number of threads.
# SNP calling done 'en masse', all processes run at same time.
# This analysis pipeline requires a lot of disk space.
# Make sure paths are set correctly as script will replace existing results directory







# Analysis functions
. functions.inc

####################### Text color variables ##############################
txtund=$(tput sgr 0 1)    # Underline
txtbld=$(tput bold)       # Bold
txtred=$(tput setaf 1)    # Red
txtgrn=$(tput setaf 2)    # Green
txtylw=$(tput setaf 3)    # Yellow
txtblu=$(tput setaf 4)    # Blue
txtpur=$(tput setaf 5)    # Purple
txtcyn=$(tput setaf 6)    # Cyan
txtwht=$(tput setaf 7)    # White
txtrst=$(tput sgr0)       # Text reset

################### Variables for analysis pipeline #########################
#------------------- Folder settings -----------------------
mosaik_dir=~andreas/mosaik-aligner/bin				#Programs bin/ directory, all programs required in same directory: Mosaik*, freebayes, bamleftalign, bamaddrg
main_dir=`pwd`								#Main directory
refseq_dir=~/galGal3/refseq						#Directory containing reference sequences
results_dir=~andreas/galGal3/results				#Directory to store the results
reads_dir=/home/gallus							#Directory containing reads
#genomes="galGal3_HS galGal3_MW"					#Name of the dataset for genome1, must correlate to BAM read-groups
genomes=`ls $reads_dir/*.dat | awk -F "." '{print $1}' | awk -F "/" '{print $4}'`			# these will be bam-file read group
# gallus_high_35_1	gallus_high_35_2	gallus_high_50	gallus_low_35_1	gallus_low_35_2	gallus_low_50

#------------------- Settings for Mosaik pipeline -----------------------
proc=16				#Number of processor cores available, -p=16
#memory=70375828		#Number of hashes to be stored in memory (default=6000000) 
#QV_filter=20		#Filtering of raw reads, phred scores, >20
hash_size=15		#Hash size, -hs=15
max_mismatches=6		#Number of allowed mismatches, -mm=4 for single-end; -mm=6 for mate-pairs
max_hash_pos=100		#Number of stored hash positions, -mhp=100
algn_thresh=25		#Alignment threshold, -act=20 for single-end, -act=25 for mate-pairs
bandwidth=17		#Bandwith for banded Smith-Waterman algorithm, -bw=13 for single-end; -bw=17 for mate-pairs
mfl=3500			#Mean fragment length in bp for mate-pairs, HARDCODED FROM ANALYSIS
radius=2000			#Search radius from mfl in bp
threads=16			#Parallel run threads in pipeline steps - Jump, Sort, DupSnoop, Merge, Coverage, Assembly

#------------------- Settings for SNP calling -----------------------
PSL=0.5			#SNP calling probability threshold
CAL=3				#SNP calling minimum number of reads
MAQ=25			#MAQ mapping alignment quality
#CRU=100			#Upper coverage limit for SNP calling
#snp_proc=34		#Number of processors available for SNP calling = maximum threads created, 34 chromosomes do all at same time


######################## Script functions ################################
Pause() {
	read -p "Press any key to continue..." -n 1 -s
}

function Presentation {
# Presentation output
echo -e "$txtgrn ########################################################"
echo -e " # A script for running Mosaik Assembler pipeline"
echo -e " # and calling SNPs and short indels using GigaBayes"
echo -e " # chromosome by chromosome for SOLiD mate pairs"
echo -e " # Author: Marcin.Kierczak@hgen.slu.se"
echo -e " # Written: 18.05.2010"
echo -e " # -------------------------------------------------"
echo -e " # Modified and expanded"
echo -e " # Co-author: Andreas.E.Lundberg@gmail.com"
echo -e " # GitHub repo: ------"
echo -e " # Date: 13.07.2011"
echo -e " ######################################################## $txtrst \n"

echo "$txtylw Program assembles one chromosome at a time from individual chromosome fasta files. \n"
echo -e " Generates individual chromosome directory in results directory. "
echo -e " First set up working directories and pipeline variables. "
echo -e " Alignment is done in linear processing fashion using multiple cores. "
echo -e " All other steps are run IN PARALLEL utilizing user selectible number of threads. "
echo -e " SNP calling done 'en masse', all processes run at same time. "
echo -e " This analysis pipeline requires a lot of disk space. "
echo -e "$txtred Make sure paths are set correctly as script will replace existing results directory !!"
echo -e " All programs need to reside in same directory !! $txtrst"
echo -e "$txtylw The reads files are processed from reads directory using list *.dat command $txtrst \n\n"
}

function DrawMenu {
echo -e "$txtblu ###########################################################################################"
echo -e " #                                       Menu                                              #"
echo -e " ########################################################################################### $txtrst \n"
echo -e "$txtcyn Initialize analysis with \t I $txtrst"
echo -e "$txtred Quit program with \t\t Q $txtrst \n"

echo -e "$txtgrn ########################### Folder paths ################################ $txtrst"
echo -e " 1) Directory of programs: \t\t\t\t\t $mosaik_dir"
echo -e " 2) Main directory for master log-file: \t\t\t $main_dir"
echo -e " 3) Reference sequence directory, one dir per chromosome: \t $refseq_dir"
echo -e " 4) Results directory: \t\t\t\t\t\t $results_dir"
echo -e " 5) Reads directory: \t\t\t\t\t\t $reads_dir"
echo -e " 6) Library name (i.e. highline & lowline): \t\t\t $genomes \n\n"


echo -e "$txtgrn ################# Assembly-specific settings ############################ $txtrst \n"
echo -e " P) Processors used for alignment, default 16: \t\t\t\t $txtred$proc$txtrst \t processors"
echo -e " H) Hash size, single-end default 15: \t\t\t\t\t $txtred$hash_size$txtrst \t hash size"
echo -e " M) Maximum mismatches allowed, defualt SE 4; defualt MP 6: \t\t $txtred$max_mismatches$txtrst \t mismatches"
echo -e " N) Maximum hash positions, default 100: \t\t\t\t $txtred$max_hash_pos$txtrst \t max hash positions"
echo -e " A) Alignment candidate threshold, SE default 20; MP default 25: \t $txtred$algn_thresh$txtrst \t alignment candidate threshold"
echo -e " B) Smith-Waterman bandwidth , SE default 13; MP default 17: \t\t $txtred$bandwidth$txtrst \t SW bandwidth"
echo -e " F) Mean fragment length/insert size: \t\t\t\t\t $txtred$mfl$txtrst\t bp mean fragment length"
echo -e " R) Search radius from mean fragmen length: \t\t\t\t $txtred$radius$txtrst\t bp radius"
echo
echo -e " T) Parallel threads created in pipeline steps, default 14: \t\t $txtred$threads$txtrst \t threads"
echo -e "    JumpDB, Sort, RemoveDuplicates, Merge, Coverage, Assembly \n\n\n"


echo -e "$txtgrn ################# Settings for SNP calling ############################## $txtrst \n\n"
echo -e " S) SNP calling probability threshold, default 0.5: \t\t\t $txtred$PSL$txtrst \t calling probabilty"
echo -e " C) Minimum number of reads per SNP call, default 3: \t\t\t $txtred$CAL$txtrst \t minimum reads"
echo -e " D) Minimum required mapping quality for alignment, default 25: \t $txtred$MAQ$txtrst \t mapping quality \n\n"

echo -e "$txtblu ################# Enter choice ############################## $txtrst \n"
echo -ne "$txtpur Enter corresponding character to update setting.$txtcyn Start analysis with I$txtred, quit with Q: $txtrst"

UpdateSettings
}

function UpdateSettings() {
	read menu_choice
	while true
	do
		case "$menu_choice" in
		"1" ) 
			echo -ne "$txtred Enter new value: $txtrst"
			read
			mosaik_dir=$REPLY
		;;
		"2" ) 
			echo -ne "$txtred Enter new value: $txtrst"
			read
			main_dir=$REPLY
		;;
		"3" ) 
			echo -ne "$txtred Enter new value: $txtrst"
			read
			refseq_dir=$REPLY
		;;
		"4" ) 
			echo -ne "$txtred Enter new value: $txtrst"
			read
			results_dir=$REPLY
		;;
		"5" ) 
			echo -ne "$txtred Enter new value: $txtrst"
			read
			reads_dir=$REPLY
		;;
		"6" ) 
			echo -ne "$txtred Enter new value: $txtrst"
			read
			genomes=$REPLY
		;;
		"P" | "p" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			proc=$REPLY
		;;
		"H" | "h" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			hash_size=$REPLY
		;;
		"M" | "m" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			max_mismatches=$REPLY
		;;
		"N" | "n" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			max_hash_pos=$REPLY
		;;
		"A" | "a" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			algn_thresh=$REPLY
		;;
		"B" | "b" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			bandwidth=$REPLY
		;;
		"F" | "f" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			mfl=$REPLY
		;;
		"R" | "r" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			radius=$REPLY
		;;
		"T" | "t" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			threads=$REPLY
		;;
		"S" | "s" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			PSL=$REPLY
		;;
		"C" | "c" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			CAL=$REPLY
		;;
		"D" | "d" )
			echo -ne "$txtred Enter new value: $txtrst"
			read
			CAL=$REPLY
		;;
		"Q" | "q" )
			echo
			exit 0
		;;
		"I" | "i" )
			break
		;;
		esac
		clear
		DrawMenu
	done
}


######################## Start script ################################

clear 
Presentation 
Pause
clear
DrawMenu

# Clean results dir, all sub-directories and log-files removed 
rm -rf $results_dir #2>/dev/null
#rm *.log 2>/dev/null

######################## Initialize pipeline #########################
echo -e "\n\n$txtpur INITIALIZE PIPELINE $txtrst"
startrun=$SECONDS
cd $refseq_dir
list=`ls $refseq_dir/v_chr*.fa | awk -F"/" '{print $NF}'`			# gives basename, i.e. chrM.fa
cd $main_dir
mkdir $results_dir 2>/dev/null

# reads=`ls $reads_dir/*.dat`



############################ Assembly loop ############################
# Re-use multi-core functions -> max threads set by $threads
for file in $list; do
	start=$SECONDS
	chromosome=${file%\.*}								# omits file extension
	mkdir $results_dir/$chromosome 2>/dev/null
	echo -e "Processing chromosome $chromosome..." | tee -a $main_dir/pipeline.log

	# Semi-sequential run
	Build	&					# detachable process
	BuildCS					# non-detachable process
	JumpDB 					# dependent on BuildCS
	wait						# wait for detached function Build to finish before entering for-loop
	
	for genome in $genomes; do
		mkdir $results_dir/$chromosome/$genome 2>/dev/null
		Align					# NO PARALLEL THREADING
	done
					
	end=$SECONDS
	exectime=$((end - start))
	echo -e "done in $exectime seconds.\n\n" | tee -a $main_dir/pipeline.log
done


#################### Alignment manipulation loop ######################
# RSMCAT, RemoveDuplicates-Sort-Merge-Coverage-Assemble-Text loop
NUM=0
QUEUE=""
for file in $list; do
	start=$SECONDS
	chromosome=${file%\.*}
	echo -e "Alignment manipulation for chromosome $chromosome..." | tee -a $main_dir/pipeline.log
	for genome in $genomes; do
		
		## All alignment manipulation loops run in parallel over $threads
		RSMCAT & 					# Start and detach process		
		PID=$!					# Get PID of process just started
		queue $PID					# 
		# Spawn process
		while [ $NUM -ge $threads ]; do 	# If $NUM is greater or equal to $threads check and regenerate queue
			checkqueue
			sleep 10
		done
	done
	end=$SECONDS
	exectime=$((end - start))
	echo -e "done in $exectime seconds.\n\n" | tee -a $main_dir/pipeline.log
done



############################# SNP calling loop ################################
# For queuing of processes
# 3 seq. runs per line -> 2 SNP calls per chromosome
NUM=0
QUEUE=""
for file in $list; do
	start=$SECONDS
	chromosome=${file%\.*}
	echo -e "Calling SNPs in chromosome $chromosome..." | tee -a $main_dir/pipeline.log
	for genome in $genomes; do
		
		## SNP calls in $threads number of threads
		CallSNP & 					# Start and detach process		
		PID=$!					# Get PID of process just started
		queue $PID					# 
		# Spawn process
		while [ $NUM -ge $threads ]; do 	# If $NUM is greater or equal to $threads check and regenerate queue
			checkqueue
			sleep 10
		done
	done
	end=$SECONDS
	exectime=$((end - start))
	echo -e "done in $exectime seconds.\n\n" | tee -a $main_dir/pipeline.log
done

endrun=$SECONDS
totaltime=$((endrun - startrun))
wait								# wait for all SNP calls to finish 

echo -e "$txtylw \n Complete assembly and SNP calling done in $totaltime seconds. \n Run completed. $txtrst \n" | tee -a $main_dir/pipeline.log




