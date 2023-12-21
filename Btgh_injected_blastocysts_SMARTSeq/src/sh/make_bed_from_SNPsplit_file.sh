#!/bin/bash

conf_threshold=${1}
SNP_file=${2}
conf_threshold_string=$(echo $conf_threshold | sed 's/\./e/')
filt_SNP_file=$(echo $SNP_file | sed "s/\.tsv/\.$conf_threshold_string\.tsv/")
genome_anno=${3}
SNP_bed=${4}
exons_SNP_bed=$(echo $SNP_bed | sed 's/\.bed/\.intersect_exons.bed/')

if grep --quiet "chr" $SNP_file
then
  awk -v conf_thr=${1} 'BEGIN {FS="\t";OFS="\t"} NR!=1 && $4>=conf_thr' $SNP_file | cat <(awk 'NR==1' $SNP_file) - > $filt_SNP_file
  awk -v conf_thr=${1} 'BEGIN {FS="\t";OFS="\t"} NR!=1 && $4>=conf_thr {print $2,$3-1,$3}' $SNP_file > $SNP_bed
else
  awk -v conf_thr=${1} 'BEGIN {FS="\t";OFS="\t"} NR!=1 && $4>=conf_thr {print $1,"chr"$2,$3,$4,$5}' $SNP_file | cat <(awk 'NR==1' $SNP_file) - > $filt_SNP_file
  awk -v conf_thr=${1} 'BEGIN {FS="\t";OFS="\t"} NR!=1 && $4>=conf_thr {print "chr"$2,$3-1,$3}' $SNP_file > $SNP_bed
fi

# intersection with exons to get gene names and to keep only SNPs overlapping exons
zcat $genome_anno | awk 'BEGIN {FS="\t";OFS="\t"} $3=="exon"' - | bedtools intersect -a $SNP_bed -b stdin -wb | awk -F";|\t" 'BEGIN {OFS="\t"} {print $1,$2,$3,$12","$15}' - | sed 's/gene_id //;s/ gene_name //' > $exons_SNP_bed
