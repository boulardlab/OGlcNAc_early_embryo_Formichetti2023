#!/bin/bash

gtf_name=$(echo $1 | sed 's/\.gtf\.gz//')

zcat $1 | awk -v FS="\t|; " '$3=="transcript" {print $10,$9}' | sed 's/"//g' | sed 's/transcript_id //' | sed 's/ gene_id /\t/' > $gtf_name.tx2gene.tsv
