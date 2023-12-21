#!/bin/bash

threads=${1}
bam_input=${2}
csorted_bam=$(echo $bam_input | sed 's/\.bam/\.csorted\.bam/' | sed 's/\.nsorted//')
bw=${3}

samtools sort -@ $threads -T tmp -o $csorted_bam $bam_input
samtools index -@ $threads $csorted_bam
bamCoverage -b $csorted_bam -o $bw -of bigwig --binSize 20 --normalizeUsing CPM -p $threads
