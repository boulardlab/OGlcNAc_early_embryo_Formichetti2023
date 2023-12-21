#!/bin/bash

threads=${1}
unassigned_bam=${2}
g1_bam=${3}
g2_bam=${4}
bam_basename=$(echo $g1_bam | sed 's/_Aligned\.out\.unique\.majchr\.nsorted\.dedupl\.genome1\.bam//')
bed_input=${5}
bedcov_output=${6}

# making new bam files after selecting only UA reads (and NOT alignment) and G1/G1 reads - respectively - from the allelic sorted alignments output by SNPsplit. This is required because for paired-end data SNPplit reports in the allelic sorted bams those alignments where at least one of the two reads is specific for one of the two genomes (and the other one can only be same genome or UA i.e. the 2 mates cannot be conflicting), while I want to use allelic-specific reads (and NOT alignments) as input for bedtools coverage 
samtools view -h $g1_bam | grep -v XX:Z:G1 | samtools view -bS - > $bam_basename.g1_UA_reads.bam
samtools view -h $g2_bam | grep -v XX:Z:G2 | samtools view -bS - > $bam_basename.g2_UA_reads.bam
samtools view -h $g1_bam | grep -v XX:Z:UA | samtools view -bS - > $bam_basename.g1_G1_reads.bam
samtools view -h $g2_bam | grep -v XX:Z:UA | samtools view -bS - > $bam_basename.g2_G2_reads.bam

# sorting by coordinates is needed for samtools merge
samtools sort -@ $threads -T tmp -o $bam_basename.g1_UA_reads.csorted.bam $bam_basename.g1_UA_reads.bam
samtools sort -@ $threads -T tmp -o $bam_basename.g2_UA_reads.csorted.bam $bam_basename.g2_UA_reads.bam
samtools sort -@ $threads -T tmp -o $bam_basename.g1_G1_reads.csorted.bam $bam_basename.g1_G1_reads.bam
samtools sort -@ $threads -T tmp -o $bam_basename.g2_G2_reads.csorted.bam $bam_basename.g2_G2_reads.bam

# indexing is necessary for samtools merge and for bedtools multicov as well
samtools index -@ $threads $bam_basename.g1_UA_reads.csorted.bam
samtools index -@ $threads $bam_basename.g2_UA_reads.csorted.bam
samtools index -@ $threads $bam_basename.g1_G1_reads.csorted.bam
samtools index -@ $threads $bam_basename.g2_G2_reads.csorted.bam
# unassigned bam from SNPsplit has been already sorted and indexed in bam_to_bw rule

# merging all bam files containing unassigned reads
echo "Merging UA reads' alignments..."
samtools merge -@ $threads $bam_basename.all_UA_reads.bam $unassigned_bam $bam_basename.g1_UA_reads.csorted.bam $bam_basename.g2_UA_reads.csorted.bam && echo "Merging finished successfully."

samtools index -@ $threads $bam_basename.all_UA_reads.bam

# computing coverage of UA, G1 and G2 reads over exon SNPs. The output will be a txt file with one column per bam input after the bed columns
#echo "Computing multicov.."
#bedtools multicov -bams $bam_basename.all_UA_reads.bam $bam_basename.g1_G1_reads.csorted.bam $bam_basename.g2_G2_reads.csorted.bam -bed $bed_input > $bedcov_output #not appropriate because it writes the coverage for the alignments and NOT the single reads!
