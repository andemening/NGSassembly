#!/bin/bash 
## Concatenate all chromosomes into whole genome single fasta file
## Covert to binary format - basespace and colorspace
## Create jump database
## by Andreas E. Lundberg 2011.06.20

mosaik_dir=~andreas/mosaik-aligner/bin
project_name=galGal3 	# project directory
refseq_dir=~/$project_name/refseq
results_dir=~/$project_name/results

## Variables
mfl=400
hash_size=15

#mkdir $refseq_dir	# dir already present containing fasta files

mkdir $results_dir 2>/dev/null
cd $refseq_dir


## Random reads moved to subdirectory refseq_random
list_random=`find $refseq_dir -name '*_random.fa'`
mkdir $refseq_dir/refseq_random 2>/dev/null
for file in $list_random; do
	mv $file $refseq_dir/refseq_random 2>/dev/null
	echo -c "Random sequences moved subdirectory refseq_random. \n $file" #| tee -a $results_dir/Jump.concat.log
done

## Start timer
start=$SECONDS

## Concat chromosomes into genome 
#Numerical order of chromosomes then alphabetical
list0=`ls chr{1..9}.fa 2>/dev/null`
list1=`ls chr1{0..9}.fa 2>/dev/null` 
list2=`ls chr2{0..9}.fa 2>/dev/null`
list3=`ls chr3{0..9}.fa 2>/dev/null`
listX=`ls chr{A..Z}*.fa 2>/dev/null`

list="$list0 $list1 $list2 $list3 $listX"
#echo $list
genome=$project_name".concat.fasta"
cat $list > $genome



## CreateBinaryFile #MosaikBuild
$mosaik_dir/MosaikBuild -fr $refseq_dir/$genome -oa $results_dir/$genome".dat" -mfl $mfl | tee -a $results_dir/Jump.concat.log

## CreateBinaryFileCS
$mosaik_dir/MosaikBuild -fr $refseq_dir/$genome -oa $results_dir/$genome".cs.dat" -mfl $mfl -cs | tee -a $results_dir/Jump.concat.log

## CreateJumpDB
st=$SECONDS
echo -e "\n\n\t\tCreating jumps database..."
# think about specifying amount of RAM used when sorting hashes with flag -mem, 25 perhaps (keys + position => ~11GB) ? NO, why? 
# flag -kd keeps keys db on disk, remove for performance
# -mhp, maximum hash positions
$mosaik_dir/MosaikJump -ia $results_dir/$genome".cs.dat" -hs $hash_size -mhp 100 -out $results_dir/$genome".cs.jmp" | tee -a $results_dir/Jump.concat.log
en=$SECONDS
exectime=$((en - st))
echo -e "\n\n\t\tJumps database done in $exectime seconds." | tee -a $results_dir/Jump.concat.log

## Fix logfile
$mosaik_dir/fixlog $results_dir/Jump.concat.log

## Calculate runtime
end=$SECONDS
exectime=$((end-start))
echo -e "\n\t\tConcat run done in $exectime seconds.\n\n" | tee -a $results_dir/Jump.concat.log

