#!/usr/bin/env bash

#SBATCH -A boulard
#SBATCH -J ggsashimi
#SBATCH -p htc-el8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH --qos=normal
#SBATCH --oversubscribe
#SBATCH -e %x-%j.out
#SBATCH -n 1
#SBATCH -o %x-%j.out

SINGULARITY="/g/boulard/Projects/embryo_public_data_reanalysis/singularity/containers/ggsashimi_latest.sif"
BAM_LIST=${1} #tsv file with col1:index of bam file, col2:path of bam file, col3:overlay level
MIN_READS_COVERAGE=${2}
GENE_NAME=${3}
#COORDINATES="chrX:101,637,000-101,686,000" #Ogt coordinates plus/minus around 1kb
COORDINATES="chr19:45,747,000-45,785,500" #Oga coordinates
OUTPUT_PATH="/g/boulard/Projects/embryo_public_data_reanalysis/analysis/sashimi/"
OUTPUT_PREFIX=$(echo $BAM_LIST | sed "s/.*\//${GENE_NAME}_/" | sed 's/.tsv//')
GTF="/g/boulard/Projects/embryo_public_data_reanalysis/data/annotations/gencode.vM25.annotation.gtf.gz"

# create output directory if it does not exist
mkdir -p $OUTPUT_PATH

# index bam files
module load SAMtools/1.16.1-GCC-11.3.0 
/g/boulard/Projects/embryo_public_data_reanalysis/src/sh/index_bam_files.sh data/sequencing/alignment/STAR/SJaccuracy/bam/sorted/

# run ggsashimi, with overlay but no aggregation
singularity run --bind /g/boulard $SINGULARITY -b $BAM_LIST -c $COORDINATES -o $OUTPUT_PATH/$OUTPUT_PREFIX -g $GTF -M $MIN_READS_COVERAGE -A median_j -O 3 -C 3 --ann-height=2 --height=2 --width=18 #rest is default
