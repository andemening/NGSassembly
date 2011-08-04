#!/bin/bash

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



