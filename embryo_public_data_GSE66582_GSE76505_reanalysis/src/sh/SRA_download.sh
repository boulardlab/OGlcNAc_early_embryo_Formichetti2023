#!/usr/bin/env bash

#SBATCH -A boulard
#SBATCH -J SRA_download
#SBATCH -p htc-el8
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH --qos=normal
#SBATCH --oversubscribe
#SBATCH -e %x-%j.out
#SBATCH -n 1
#SBATCH -o %x-%j.out

SRA_ACC_LIST_FILE=$1
SEQ_DATA_DIR=$2
SRA_DIR_PREFIX=$SEQ_DATA_DIR/sra
FASTQ_DIR=$SEQ_DATA_DIR/fastq/
TMP_DIR=$SRA_DIR_PREFIX/tmp
#CONDA_ENV=/g/boulard/sformich/conda/envs/smake_SRA_download

# create output directories if they do not exist
mkdir -p $SRA_DIR_PREFIX
mkdir -p $FASTQ_DIR

while IFS=$'\n' read -r SRA_Acc; do
    wget --directory-prefix=$SRA_DIR_PREFIX/ https://sra-pub-run-odp.s3.amazonaws.com/sra/"$SRA_Acc"/"$SRA_Acc"
    mkdir -p $TMP_DIR/fastqdump-"$SRA_Acc"
  parallel-fastq-dump -s $SRA_DIR_PREFIX/"$SRA_Acc" -t 20 -O $FASTQ_DIR --tmpdir $TMP_DIR/fastqdump-"$SRA_Acc" --split-files
    gzip "$SRA_Acc"_1.fastq
    gzip "$SRA_Acc"_2.fastq 
    rm -r $TMP_DIR/fastqdump-"$SRA_Acc"
done < $SRA_ACC_LIST_FILE
