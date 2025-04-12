# Formichetti et al., 2024 - Perturbing nuclear glycosylation in the mouse preimplantation embryo slows down embryonic growth

# Introduction

* This repository contains all code used to generate the figures and conclusions based on RNA-Seq data contained in paper "Perturbing nuclear glycosylation in the mouse preimplantation embryo slows down embryonic growth" ([Formichetti et al. 2025](https://www.pnas.org/doi/10.1073/pnas.2410520122)).

* For each sequencing dataset, there is a subrepository with a self-explanatory name: 

    - 4 repositories for the 4 Smart-Seq datasets generated in our study: 2-cell embryos, morulae, blastocysts and E7 embryos non-injected or injected with Btgh/dBtgh
    - 1 repository for all analyses performed using publicly available RNA-Seq data spanning mouse embryonic development (from [GSE66582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66582) and [GSE76505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76505)).

# General overview of the workflows

## 2-cell embryos, morulae, blastocysts and E7 embryos from this study

For each of the 4 embryonic stages, the general workflow is the same:

1. a pipeline starting from fastq files to read counts, which was:

    - for single copy genes: a Galaxy pipeline, found in **src/galaxy** of each subrepository
    - for retrotrasposons: a custom snakemake pipeline, found in **snake-make/TE_RNASeq.Snakefile** of each subrepository
    - for the allele-specific RNA-Seq analysis of blastocysts: a custom snakemake pipeline, found in **snake-make/SNPsplit.Snakefile**
    - for the snakemake pipelines, config files are in config/ and conda environments are in env/conda

<n>

2. different kinds of custom downstream analyses using the output of each pipeline and performed with R, all included in Rmd files with self-explanatory names, found in **src/Rmd** and **whose output is [here](https://boulardlab.github.io/OGlcNAc_early_embryo_Formichetti2023/)**. Order of Rmd is specified in the main README of each subdirectory.

For **2-cell embryos**, the repository also contains the analysis of Smart-Seq data from [GSE111864](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSE111864) (the subset of samples analyzed is in sequencing/Conine/SRR_Acc_List_cauda_2C_GSE111864.txt):

- same Galaxy pipeline used for the analysis of single copy genes in the Smart-Seq data generated in our study
- downstream analysis in a custom Rmd (src/Rmd/Conine.Rmd) whose outputs (list of maternal and EGA genes and counts table) are used in src/Rmd/embryo_filtering_and_clustering.Rmd

## GSE66582 and GSE76505

### Analysis of expression:

1. a pipeline starting from fastq files to read counts, which was:

    - for single copy genes: a Galaxy pipeline, found in **src/galaxy/Galaxy-Workflow-PE_TranscriptCount_v1.3.ga** (which outputs transcripts counts, summarized at gene level downstream in Rmd)
    - for retrotransposons: a custom snakemake pipeline, found in **snake-make/TE_RNASeq.Snakefile**

<n>

2. different kinds of custom downstream analyses using the output of each pipeline and performed with R, all included in Rmd files with self-explanatory names, found in **src/Rmd** and **whose output is [here](https://boulardlab.github.io/OGlcNAc_early_embryo_Formichetti2023/)**. Order of Rmd is specified in the main README of the subdirectory.

### Analysis of reads coverage / alternative splicing:

1. a Galaxy pipeline starting from fastq files and producing the bigwig files. It uses more sensitive STAR mapping parameters than the one used for analysis of expression and can be found in **src/galaxy/Galaxy-Workflow-PE_SMART-Seq_embryo_mapToGenome_SJaccuracy.ga**  

2. to produce sashimi plots for *[Ogt](https://github.com/boulardlab/OGlcNAc_early_embryo_Formichetti2023/blob/main/embryo_public_data_GSE66582_GSE76505_reanalysis/analysis/sashimi/pdf/Ogt_bam_list.pdf)* and *[Oga](https://github.com/boulardlab/OGlcNAc_early_embryo_Formichetti2023/blob/main/embryo_public_data_GSE66582_GSE76505_reanalysis/analysis/sashimi/pdf/Oga_bam_list.pdf)*, [ggsashimi](https://github.com/guigolab/ggsashimi) (called in src/sh/ggsashimi.sh).

# Raw data

Sequencing data generated in our study are available at [E-MTAB-12981](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12981?key=29a1ff26-6018-4d0c-a08d-21d23226cb65), which contains 4 sequencing run, one for each of the 4 embryonic stages analyzed after injection of Btgh/dBtgh/non-injected.






