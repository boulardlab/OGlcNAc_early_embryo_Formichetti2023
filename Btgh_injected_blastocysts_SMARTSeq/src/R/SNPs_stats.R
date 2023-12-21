#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# creating output directory if it does not exist
system("[ -d ../../analysis/SNPs_stats/ ] || mkdir -p ../../analysis/SNPs_stats/")

# load transcript lengths derived from gencode GTF
tx_len <- fread(snakemake@config[["tx_lengths"]])
names(tx_len) <- c("gene_id", "gene_name", "tx_length")
# summarizing at gene level by taking length of longer isoform for each gene as well as minimum and avg length of isoforms
tx_len[ , c("max_len", "min_len", "avg_len") := list(max(tx_length), min(tx_length), mean(tx_length)), by=.(gene_id, gene_name)]
# removing version from gene id
tx_len$gene_id <- gsub("\\..*", "", tx_len$gene_id)

# loading bed file of exons-only SNPs (with gene info) obtained using chosen confidence threshold
bed_exonsOnly <- fread(snakemake@input[["exons_SNP_bed"]], sep = "\t")
bed_exonsOnly$gene_id <- gsub("\\..*", "", bed_exonsOnly$V4)
bed_exonsOnly$gene_name <- gsub('.*,"', "", bed_exonsOnly$V4)
bed_exonsOnly$V4 <- NULL
names(bed_exonsOnly) <- c("chr", "start", "end", "gene_id", "gene_name")
# loading tsv file containing all SNPs and confidence level for PWD of each of them
tsv <- fread(snakemake@config[["SNP_tsv"]])
# merging it with the bed in order to look only at SNPs in exons
bed_exonsOnly$Chr <- gsub("chr", "", bed_exonsOnly$chr)
tsv_exonsOnly <- merge(tsv, bed_exonsOnly, by.x = c("Chr", "Position"), by.y = c("Chr", "end"))

## plotting histogram of confidence levels
pdf(snakemake@output[["conf_level_hist"]])
ggplot(data = tsv, aes(x = `SNP value`)) +
  geom_histogram(bins = 40) +
  ggtitle("all FVB-confident SNPs")
ggplot(data = tsv_exonsOnly, aes(x = `SNP value`)) +
  geom_histogram(bins = 40) +
  ggtitle("only FVB-confident SNPs in exons")
dev.off()

## function which - from the bed file containing also the gene info - calculates the number of SNPs per gene and adds the gene lengths info. Then it summarizes the table at gene level removing SNPs columns, for plotting and stats purposes
make_gene_stats_dt <- function (my_bed) {
  # removing chrM genes
  my_bed <- my_bed[chr!="chrM"]
  # computing number of SNPs per gene
  my_bed[, SNPs_per_gene := .N, by=.(gene_id, gene_name, Chr)]
  # adding tx lengths info
  my_bed <- merge(my_bed, unique(tx_len[,.(gene_id, max_len, min_len, avg_len)]), by="gene_id")
  # summarizing table at gene level
  genes_stats_dt <- unique(my_bed[,.(gene_id, chr, gene_name, SNPs_per_gene, max_len, min_len, avg_len)])
  # computing SNPs per kb, using both longer and shorter RNA isoform for each gene 
  genes_stats_dt$max_SNPs_per_kb <- genes_stats_dt$SNPs_per_gene / genes_stats_dt$min_len * 1000
  genes_stats_dt$min_SNPs_per_kb <- genes_stats_dt$SNPs_per_gene / genes_stats_dt$max_len * 1000
  return(genes_stats_dt)
}

gene_SNP_stats <- make_gene_stats_dt(bed_exonsOnly)

## plotting SNPs per gene
pdf(snakemake@output[["SNPs_per_gene_hist"]])
ggplot(data = gene_SNP_stats, aes(x = SNPs_per_gene)) +
  geom_histogram(binwidth = 1, boundary = 0) + 
  scale_x_continuous(breaks=c(0, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50), limits=c(-1,51)) +
  ggtitle("conf threshold = 0")
dev.off()

## plotting SNPs per kb
pdf(snakemake@output[["SNPs_per_kb_hist"]])
ggplot(data = gene_SNP_stats, aes(x = max_SNPs_per_kb)) +
  geom_histogram(binwidth = 1, boundary = 0) + 
  scale_x_continuous(breaks=c(0, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50), limits=c(-1,51)) +
  ggtitle("taking min tx isoform length")
ggplot(data = gene_SNP_stats, aes(x = min_SNPs_per_kb)) +
  geom_histogram(binwidth = 1, boundary = 0) + 
  scale_x_continuous(breaks=c(0, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50), limits=c(-1,51)) +
  ggtitle("taking max tx isoform length")
dev.off()
