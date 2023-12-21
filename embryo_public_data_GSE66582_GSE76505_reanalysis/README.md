Analysis of public datasets from [GSE66582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66582) (MII oocyte to inner cell mass) and [GSE76505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76505) (from E3.5 to E7.5). Only WT samples were analyzed and they are listed in data/sequencing/SraRunTable_selection.csv (for the analysis of single copy genes) and data/sequencing/SraRunTable_selection_TE_analysis.csv. 

# Order of Rmd

## Single copy genes

All Rmd start from the output of galaxy pipeline in src/galaxy/Galaxy-Workflow-PE_TranscriptCount_v1.3.ga.

## Retrotransposons

TE_analysis.Rmd starts from the output of TE_RNASeq.Snakefile.
