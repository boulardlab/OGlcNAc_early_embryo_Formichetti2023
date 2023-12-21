#!/usr/bin/env bash

for INFILE in $(ls $1/*.bam)
do
   samtools index $INFILE -@ 2
done
