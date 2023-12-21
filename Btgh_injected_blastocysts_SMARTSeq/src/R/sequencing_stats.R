#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

seq_stat_files <- c("STAR_input_reads.txt", "featureCounts_assigned.txt", "STAR_uniquely_mapped.txt", "STAR_multiple_loci.txt", "STAR_too_many_loci.txt")
names(seq_stat_files) <- gsub("\\.txt", "", seq_stat_files)
STAR_default_folder <- "../../data/sequencing/stats/STAR_default/"
STAR_TE_folder <- "../../data/sequencing/stats/STAR_TE/"
STAR_riboDepl_TE_noIntrons_folder <- "../../data/sequencing/stats/STAR_riboDepl_TE/noIntrons/"
STAR_riboDepl_TE_wIntrons_folder <- "../../data/sequencing/stats/STAR_riboDepl_TE/wIntrons/"

get_stat_list <- function (my_stat_folder, i, my_stat_files) {
  # by using fread w these sep and sep2 parameters, the numeric column is always V2
  stat_dt <- fread(paste0(my_stat_folder, my_stat_files[i]), sep = "\t", sep2 = " ")
  names(stat_dt) <- c("sample", names(my_stat_files)[i])
  stat_dt$sample <- gsub(".*lane1", "", gsub("\\.txt.*", "", stat_dt$sample))
  stat_dt$sample <- gsub("_.*", "", gsub(".*/", "", stat_dt$sample))
  return(stat_dt)
}

make_stat_plots <- function (stat_folder, stat_files = seq_stat_files) {
  
  stat_list <- lapply(1:length(stat_files), get_stat_list, my_stat_folder = stat_folder, my_stat_files = stat_files)
  
  all_stat_dt <- Reduce(function(...) merge(..., by = c("sample")), stat_list)
  my_cols <- c("sample", names(stat_files)[-1])
  all_stat_melt_dt <- melt(all_stat_dt[, ..my_cols], variable.name = "stat")
  my_cols <- c("sample", "STAR_input_reads")
  all_stat_melt_dt<- merge(all_stat_melt_dt, all_stat_dt[, ..my_cols], by = "sample")
  
  pdf(file = paste0(stat_folder, "stats.pdf"))
  num_violin <- ggplot(data = all_stat_melt_dt, aes(x = stat, y = value)) +
                  geom_violin() +
                  stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black") +
                  ylab("number")
  perc_violin <- ggplot(data = all_stat_melt_dt, aes(x = stat, y = value/STAR_input_reads*100)) +
                  geom_violin() +
                  stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black") +
                  ylab("perc")
  print(num_violin)
  print(perc_violin)
  dev.off()
  
}

#make_stat_plots(stat_folder = STAR_default_folder)
make_stat_plots(stat_folder = STAR_TE_folder)
#make_stat_plots(stat_folder = STAR_riboDepl_TE_noIntrons_folder)
#make_stat_plots(stat_folder = STAR_riboDepl_TE_wIntrons_folder)


