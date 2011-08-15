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
	$mosaik_dir/bamaddrg $infiles | $mosaik_dir/freebayes --stdin -f $refseq_dir/$file --pooled --ploidy $ploidy -dd --pvar $poly_prob --no-filters --min-coverage $CAL --min-mapping-quality $MAQ --ignore-reference-allele --use-mapping-quality -v $results_dir/$file.called.vcf

# Merge bam files in stream and pipe to BEDtools genomeCoverageBed for coverage calculation of histogram  
~/bamtools/bin/bamtools merge -in ~/NGSassembly/gallus.results/SNPcall/gallus_high_frag_35bp_1.filtered.dat.aligned.sorted.Marked.bam -in ~/NGSassembly/gallus.results/SNPcall/gallus_high_frag_35bp_2.filtered.dat.aligned.sorted.Marked.bam -in ~/NGSassembly/gallus.results/SNPcall/gallus_high_mate_2x50bp.filtered.dat.aligned.sorted.Marked.bam | ~/BEDtools/bin/genomeCoverageBed -ibam stdin -bga -g chromInfo.txt > coverage.galGal3.high.bedgraph

Output stream from BEDtools: bamtools merge | BEDtools genomeCoverageBed | awk calculation, no saving of histogram coverage file
~/bamtools/bin/bamtools merge -in ~/NGSassembly/gallus.results/SNPcall/gallus_high_frag_35bp_1.filtered.dat.aligned.sorted.Marked.bam -in ~/NGSassembly/gallus.results/SNPcall/gallus_high_frag_35bp_2.filtered.dat.aligned.sorted.Marked.bam -in ~/NGSassembly/gallus.results/SNPcall/gallus_high_mate_2x50bp.filtered.dat.aligned.sorted.Marked.bam | ~/BEDtools/bin/genomeCoverageBed -ibam stdin -d -g ../chromInfo.no_random.txt 


# Calculation of X coverage using samtools merge piped to BEDtools genomeCoverageBed histogram file (all bases covered)
# Combined total bases, mean coverage and max
cat filename | awk -F "\t" '{print $3}' | awk '{for (i=1; i<=NF; i++) s=s+$i} $1>=max {max=$1} END {print "\nTotal bases: "NR, "\nMean coverage: "s/NR, "\nMaximum coverage (minimum coverage is 0): "max}'

# with SD 
time awk -F "\t" '{print $3}' filename | awk '{for (i=1; i<=NF; i++) s=s+$i} {sum+=$1;sumsq+=$1*$1} $1>=max {max=$1} END {print "\nTotal bases: "NR, "\nMean coverage: "s/NR, "(SD: "(sqrt(sumsq/NR-(sum/NR)^2)), ")", "\nMaximum coverage (minimum coverage 0): "max}'

# With % coverage, how many percent of genome is covered with contigs
time awk -F "\t" '{print $3}' filename | awk '$1!=0 {contigs++} {for (i=1; i<=NF; i++) s=s+$i} {sum+=$1;sumsq+=$1*$1} $1>=max {max=$1} END {print "\n\nTotal bases: "NR, "\nMean coverage: "s/NR, "( SD: "(sqrt(sumsq/NR-(sum/NR)^2)), ")", "\nMaximum coverage: "max, "( minimum coverage 0 )", "\n\nBases with coverage / contigs: "contigs, "\nContig coverage, percent of genome covered by contigs: "contigs*100/NR, "%"}'

Everything at the same time: 
~/bamtools/bin/bamtools merge -in ~/NGSassembly/gallus.results/SNPcall/gallus_high_frag_35bp_1.filtered.dat.aligned.sorted.Marked.bam -in ~/NGSassembly/gallus.results/SNPcall/gallus_high_frag_35bp_2.filtered.dat.aligned.sorted.Marked.bam -in ~/NGSassembly/gallus.results/SNPcall/gallus_high_mate_2x50bp.filtered.dat.aligned.sorted.Marked.bam \
| ~/BEDtools/bin/genomeCoverageBed -ibam stdin -d -g ../chromInfo.no_random.txt \
|& tee -a bajs | time awk -F "\t" '{print $3}' \
| awk '$1!=0 {contigs++} {for (i=1; i<=NF; i++) s=s+$i} {sum+=$1;sumsq+=$1*$1} $1>=max {max=$1} END {print "\n\nTotal bases: "NR, "\nMean coverage: "s/NR, "( SD: "(sqrt(sumsq/NR-(sum/NR)^2)), ")", "\nMaximum coverage: "max, "( minimum coverage 0 )", "\n\nBases with coverage / contigs: "contigs, "\nContig coverage, percent of genome covered by contigs: "contigs*100/NR, "%"}'