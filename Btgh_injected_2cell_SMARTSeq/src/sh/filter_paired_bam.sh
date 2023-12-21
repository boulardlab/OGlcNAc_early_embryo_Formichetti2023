#!/bin/bash

threads=${1}
TE_annotation=${2}
gene_annotation=${3}
bam_input=${4}
tmp_dir=${5}

sample_name=$(echo $bam_input | sed 's/_Aligned.*//' | sed 's/.*\///')
bam_dir=$(echo $bam_input | sed 's/\/[a-zA-Z0-9_.]*$//')

# keeping only reads aligning inside TE; in fact, complete overlap is required by the -f parameter. Of these, excluding reads overlapping with gene bodies
bedtools intersect -abam $bam_input -b $TE_annotation -u -f 1 | bedtools intersect -abam stdin -b $gene_annotation -v | samtools view -@ $threads -o $tmp_dir$sample_name.sam

# finding those reads for which only one of the two mates passed the filters in above line and writing their qnames
cut -f1 $tmp_dir$sample_name.sam | sort | uniq -c | awk '$1!=2{print $2}' > $tmp_dir$sample_name.qnames.txt

# removing these reads from the sam file, reattaching SAM header and converting back to BAM format 
grep -F -v -f $tmp_dir$sample_name.qnames.txt $tmp_dir$sample_name.sam | cat <(samtools view -H $bam_input) - | samtools view -@ $threads -b - -o $bam_dir/$sample_name.TE_only.bam

#wc -l $tmp_dir$sample_name.sam > counts.txt
#wc -l $tmp_dir$sample_name.nof.sam >> counts.txt
#wc -l $tmp_dir$sample_name.qnames.txt >> counts.txt
#wc -l $tmp_dir$sample_name.qnames.nof.txt >> counts.txt
#samtools view $bam_dir/$sample_name.TE_only.bam | wc -l >> counts.txt
#samtools view $bam_dir/$sample_name.TE_only.nof.bam | wc -l >> counts.txt
