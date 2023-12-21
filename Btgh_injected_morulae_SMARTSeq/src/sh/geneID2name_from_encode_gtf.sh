#!/bin/bash

gtf_name=$(echo $1 | sed 's/\.gtf\.gz//')

zcat $1 | awk -v FS="\t|; " '$3=="gene" {print $9,$11}' | sed 's/"//g' | sed 's/gene_id //' | sed 's/ gene_name /\t/' > $gtf_name.geneID2name.tsv
