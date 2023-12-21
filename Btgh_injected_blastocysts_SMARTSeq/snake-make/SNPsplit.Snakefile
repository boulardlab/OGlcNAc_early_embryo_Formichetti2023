########################################################################################
# Pipeline for analysis of differential allelic expression in RNA-Seq of F1 hybrid samples
# Date: 17/08/22
# Author: Sara Formichetti
# email: sara.formichetti@embl.it
# run: all specified in shell wrapper for the specific experiment's analysis
#######################################################################################

#####################
# Imports
#####################

import pandas
import os
import yaml

#####################
# Defining shell to use
#####################

shell.executable("/bin/bash")

#####################
# Create tmp dir if it does not exist
#####################
if not os.path.isdir('./tmp'):
    shell("mkdir tmp/")

#####################
# Functions
#####################

def read_samplesTable(samplesTable):
    samplesData = pandas.read_csv(samplesTable)
    # Verify column names
    if not {'SRR','sample','condition', 'group_or_time_point', 'batch', 'biol_rep','tech_rep','run','library_layout','read_length'}.issubset(samplesData.columns.values):
	    raise KeyError("The samples file must contain the following named columns (name of the column[content format]): SRR[SRR ID i.e. string],sample[string with no spaces],condition[string with no spaces],group_or_time_point[string with no spaces],batch[number],biol_rep[number],tech_rep[number],run[number],library_layout[SINGLE - |PAIRED - ],read_length[number]; please add this column and just write always 'none'(without quotes) under 'SRR' if they do not come from SRA and always '1'(without quotes) under 'condition/group_or_time_point/biol_rep/tech_rep' if there is no more than one group.")
    return samplesData

######################
# Config Variables
######################

# Checking that provided configfile contains all necessary parameters 

with open('config/SNPsplit_config.yaml') as f:
    my_config_dict = yaml.safe_load(f)

target_config_keys = ['experiment_folder', 'input_samples_table', 'fastq_input_dir', 'seq_output_dir', 'genome_dir', 'genome_fasta', 'SNP_tsv', 'conf_threshold', 'maj_chr_bed', 'annotation_dir', 'gene_annotation_bed', 'gene_annotation_gtf', 'tx_lengths', 'analysis_output_dir', 'logs_dir', 'tmp_dir', 'envs_folder', 'script_dir']
missing_params = [ele for ele in target_config_keys if ele not in list(my_config_dict)] 
if len(missing_params) == 0:
    print("Configfile contains all necessary variables")
else :
    print("Configfile misses following parameters: ", str(missing_params))
    sys.exit()

# Getting value of variable containing samples' table file name

try:
  input_samples_table = config["input_samples_table"]
except KeyError:
    print("The parameter \"input_samples_table\" has not been defined in the configfile. The pipeline cannot start.")
    
#################################
# Reading input samples' table
#################################

inputSamplesData = read_samplesTable(input_samples_table)

########################
# Variables definition
########################

# Splitting the table into single or paired end experiments

index_single = inputSamplesData['library_layout'] == 'SINGLE - '
index_paired = inputSamplesData['library_layout'] == 'PAIRED - '
samplesData_single = inputSamplesData[index_single]
samplesData_paired = inputSamplesData[index_paired]

# Output files names

SINGLESAMPLES = samplesData_single['sample'].tolist()
PAIREDSAMPLES = samplesData_paired['sample'].tolist()

###############################################################################
# Rules
###############################################################################

######################
# Rule all
######################

rule all:
  input:
   #expand("{path}qc/fastqc/{sample}_{mate}.fastqc.{ext}", path=config["seq_output_dir"], sample=PAIREDSAMPLES, ext=["html", "zip"], mate=["1","2"]),
    os.path.join(config["analysis_output_dir"], "SNPs_stats/SNPs_per_kb.pdf"),
   #expand("{path}alignment/STAR/maskedGenome/SNPsplit/bedcov/{sample}_Aligned.out.unique.majchr.dedupl.split.SNP_cov.txt", path=config["seq_output_dir"], sample=PAIREDSAMPLES),
    expand("{path}featureCounts/SNPsplit/{sample}_joined_counts.txt", path=config["seq_output_dir"], sample=PAIREDSAMPLES),
    os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/STAR_multiple_loci.txt"),
    os.path.join(config["seq_output_dir"], "stats/SNPsplit/g2_reads.txt"),
    os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/featureCounts_assigned_g2.txt")

rule fastqc_my_fastq:
  input:
    expand("{path}{{sample}}_{{mate}}.fastq.gz", path=config["fastq_input_dir"])
  params:
    output_dir=expand("{path}qc/fastqc/", path=config["seq_output_dir"])
  output:
    expand("{path}qc/fastqc/{{sample}}_{{mate}}.fastqc.{ext}", path=config["seq_output_dir"], ext=["html", "zip"])
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml" 
  shell:
    "fastqc -t 6 -o {params.output_dir} {input}"

rule sort_my_fastq:
  input:
    expand("{path}{{sample}}_{{mate}}.fastq.gz", path=config["fastq_input_dir"])
  output:
    expand("{path}{{sample}}_{{mate}}.sorted.fastq.gz", path=config["fastq_input_dir"])
  threads: 6
  shell:
    """
    zcat {input} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | gzip - > {output}
    """

rule trimgalore_my_fastq:
  input:
    R1=expand("{path}{{sample}}_1.sorted.fastq.gz", path=config["fastq_input_dir"]),
    R2=expand("{path}{{sample}}_2.sorted.fastq.gz", path=config["fastq_input_dir"])
  output:
    R1_trimmed_fq=expand("{path}trimming/{{sample}}_1.sorted_val_1.fq.gz", path=config["seq_output_dir"]),
    R2_trimmed_fq=expand("{path}trimming/{{sample}}_2.sorted_val_2.fq.gz", path=config["seq_output_dir"])
  params:
    output_dir=expand("{path}trimming/", path=config["seq_output_dir"]),
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell: 
    "trim_galore -j 6 -q 20 --stringency 1 -e 0.1 --length 20 --fastqc -o {params.output_dir} --paired {input.R1} {input.R2}" # trimgalore default parameters

rule get_SNP_bed:
  input:
    SNP_file=config["SNP_tsv"],
    genome_annotation=config["gene_annotation_gtf"]
  output:
    filtered_SNP_file=config["SNP_tsv"].replace(".tsv", ".") + str(config["conf_threshold"]).replace(".", "e") + ".tsv",
    SNP_bed=config["SNP_tsv"].replace("_w_header.tsv", ".") + str(config["conf_threshold"]).replace(".", "e") + ".bed",
    exons_SNP_bed=config["SNP_tsv"].replace("_w_header.tsv", ".") + str(config["conf_threshold"]).replace(".", "e") + ".intersect_exons.bed"
  params:
    conf_threshold=config["conf_threshold"]
  threads: 4
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "{config[script_dir]}sh/make_bed_from_SNPsplit_file.sh {params.conf_threshold} {input.SNP_file} {input.genome_annotation} {output.SNP_bed}" 

rule SNPs_stats:
  input:
    exons_SNP_bed=rules.get_SNP_bed.output.exons_SNP_bed
  output:
    conf_level_hist=os.path.join(config["analysis_output_dir"], "SNPs_stats/exons_SNPs_conf_levels.pdf"),
    SNPs_per_gene_hist=os.path.join(config["analysis_output_dir"], "SNPs_stats/SNPs_per_gene.pdf"),
    SNPs_per_kb_hist=os.path.join(config["analysis_output_dir"], "SNPs_stats/SNPs_per_kb.pdf")
  threads: 4
  conda: config["envs_folder"] + "SNPsplit_R.yml"
  script:
    config["experiment_folder"] + config["script_dir"] + "R/SNPs_stats.R"

rule mask_genome:
  input:
    genome_fa=config["genome_fasta"],
    SNP_bed=rules.get_SNP_bed.output.SNP_bed
  output:
    tmp_unz_fa=temp(config["genome_fasta"].replace(".fa.gz", ".fa")),
    unz_masked_fa=temp(config["genome_fasta"].replace(".fa.gz", ".masked.fa")),
    masked_fa=config["genome_fasta"].replace(".fa.gz", ".masked.fa.gz")
  threads: 4
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    """
    zcat {input.genome_fa} > {output.tmp_unz_fa}
    bedtools maskfasta -fi {output.tmp_unz_fa} -bed {input.SNP_bed} -fo {output.unz_masked_fa}
    gzip -c {output.unz_masked_fa} > {output.masked_fa}
    """

rule STAR_genome_preparation:
  input:
    unz_masked_fa=rules.mask_genome.output.unz_masked_fa,
    genome_annotation=config["gene_annotation_gtf"]
  output:
    unz_gtf=temp(config["gene_annotation_gtf"].replace(".gtf.gz", ".gtf")),
    index_dir=directory(os.path.join(config["genome_dir"], "STAR/index_2.7.5c"))
  params:
    tmp_folder=os.path.join(config["tmp_dir"], "STAR_genome_prep")
  log:
    os.path.join(config["logs_dir"], "STAR_genome_prep.log")	  
  threads: 8
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    """
    ## tmp directory needs to be removed before running genomeGenerate otherwise STAR complains
    [[ -d {params.tmp_folder} ]] && rm -rf {params.tmp_folder}
    ## STAR wants uncompressed genome fasta
    zcat {input.genome_annotation} > {output.unz_gtf}
    STAR --runMode genomeGenerate \
    --outTmpDir {params.tmp_folder} \
    --runThreadN {threads} \
    --genomeDir {output.index_dir} \
    --genomeFastaFiles {input.unz_masked_fa} \
    --sjdbGTFfile {output.unz_gtf} \
    --sjdbOverhang 39 |& \
    tee {log}
    """

rule STAR_align:
  input:
    R1=rules.trimgalore_my_fastq.output.R1_trimmed_fq,
    R2=rules.trimgalore_my_fastq.output.R2_trimmed_fq,
    STAR_genome_dir=rules.STAR_genome_preparation.output.index_dir
  output:
    bam=expand("{path}alignment/STAR/maskedGenome/{{sample}}_Aligned.out.bam", path=config["seq_output_dir"])
  params:
    output_prefix=expand("{path}alignment/STAR/maskedGenome/{{sample}}_", path=config["seq_output_dir"])
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "STAR --runThreadN {threads} --genomeDir {input.STAR_genome_dir} --readFilesCommand zcat --readFilesIn {input.R1} {input.R2} --outFileNamePrefix {params.output_prefix} --runMode alignReads --outSAMtype BAM Unsorted --alignEndsType EndToEnd --outSAMattributes NH HI NM MD --outSAMprimaryFlag OneBestScore --outSAMmapqUnique 255 --outSAMunmapped None --outFilterIntronStrands RemoveInconsistentStrands --outFilterIntronMotifs RemoveNoncanonical --outFilterType Normal --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 1.0 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 --outSAMmultNmax -1 --outSAMtlen 1" ### except for --alignEndsType and --outSANattributes values - required by SNPsplit - the rest is default
## the parameters controlling allowed mismatches (mainly --outFilterMismatchNmax) are fine. In fact, in script SNPs_stats.R I computed distribution of number of SNPs per gene per kb and the 75% percentile is around 13 i.e. from 1 to 2 SNPs per read pair, which is way lower 10. 

rule get_uniquely_al_bam:
  input:
    rules.STAR_align.output.bam
  output:
    expand("{path}alignment/STAR/maskedGenome/filtered/{{sample}}_Aligned.out.unique.bam", path=config["seq_output_dir"])
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "samtools view -b -q 255 {input} > {output}"
## using only uniquely aligned reads is suggested in Castel et al., Genome Biol 2015

rule nsort_maj_chr_bam:
  input:
    rules.get_uniquely_al_bam.output
  output:
    expand("{path}alignment/STAR/maskedGenome/filtered/{{sample}}_Aligned.out.unique.majchr.nsorted.bam", path=config["seq_output_dir"])
  params:
    majchr_bed=config["maj_chr_bed"]
  threads: 6 
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "samtools view -@ {threads} -L {params.majchr_bed} {input} -b | samtools sort -n -@ {threads} -T tmp -o {output}"

rule remove_duplicates:
  input:
    rules.nsort_maj_chr_bam.output
  output:
    dupl_rm_bam=expand("{path}alignment/STAR/maskedGenome/filtered/{{sample}}_Aligned.out.unique.majchr.nsorted.dedupl.bam", path=config["seq_output_dir"]),
    dupl_metrics=expand("{path}alignment/STAR/maskedGenome/filtered/{{sample}}_dedupl_metrics.txt", path=config["seq_output_dir"])
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "picard MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=queryname I={input} O={output.dupl_rm_bam} M={output.dupl_metrics}" # by default, it chooses the read to keep with an algorithm that ranks reads by the sums of their base-quality scores
## duplicates removal is suggested in Castel et al., Genome Biol 2015

rule SNPsplit_bam:
  input: 
    bam=rules.remove_duplicates.output.dupl_rm_bam,
    SNP_file=rules.get_SNP_bed.output.filtered_SNP_file
  output:
    flagged_bam=expand("{path}alignment/STAR/maskedGenome/SNPsplit/{{sample}}_Aligned.out.unique.majchr.nsorted.dedupl.allele_flagged.bam", path=config["seq_output_dir"]),
    unassigned_bam=expand("{path}alignment/STAR/maskedGenome/SNPsplit/{{sample}}_Aligned.out.unique.majchr.nsorted.dedupl.unassigned.bam", path=config["seq_output_dir"]),
    g1_bam=expand("{path}alignment/STAR/maskedGenome/SNPsplit/{{sample}}_Aligned.out.unique.majchr.nsorted.dedupl.genome1.bam", path=config["seq_output_dir"]),
    g2_bam=expand("{path}alignment/STAR/maskedGenome/SNPsplit/{{sample}}_Aligned.out.unique.majchr.nsorted.dedupl.genome2.bam", path=config["seq_output_dir"])
  params:
    outdir=os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/")
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "SNPsplit --paired -o {params.outdir} --snp_file {input.SNP_file} {input.bam}"

rule flagged_bam_to_bw:
  input:
    rules.SNPsplit_bam.output.flagged_bam
  output:
    csorted_bam=os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/{sample}_Aligned.out.unique.majchr.dedupl.allele_flagged.csorted.bam"),
    bw=os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/bw/{sample}_Aligned.out.unique.majchr.dedupl.allele_flagged.bw")
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "{config[script_dir]}sh/bam_to_bw.sh {threads} {input} {output.bw}"

rule unassigned_bam_to_bw:
  input:
    rules.SNPsplit_bam.output.unassigned_bam
  output:
    csorted_bam=os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/{sample}_Aligned.out.unique.majchr.dedupl.unassigned.csorted.bam"),
    bw=os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/bw/{sample}_Aligned.out.unique.majchr.nsorted.dedupl.unassigned.bw")
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "{config[script_dir]}sh/bam_to_bw.sh {threads} {input} {output.bw}"

rule g1_bam_to_bw:
  input:
    rules.SNPsplit_bam.output.g1_bam
  output:
    csorted_bam=os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/{sample}_Aligned.out.unique.majchr.dedupl.genome1.csorted.bam"),
    bw=os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/bw/{sample}_Aligned.out.unique.majchr.nsorted.dedupl.genome1.bw")
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "{config[script_dir]}sh/bam_to_bw.sh {threads} {input} {output.bw}"

rule g2_bam_to_bw:
  input:
    rules.SNPsplit_bam.output.g2_bam
  output:
    csorted_bam=os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/{sample}_Aligned.out.unique.majchr.dedupl.genome2.csorted.bam"),
    bw=os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/bw/{sample}_Aligned.out.unique.majchr.nsorted.dedupl.genome2.bw")
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "{config[script_dir]}sh/bam_to_bw.sh {threads} {input} {output.bw}"

rule split_reads_SNP_coverage:
  input:
    unassigned_bam=rules.unassigned_bam_to_bw.output.csorted_bam, 
    g1_bam=rules.SNPsplit_bam.output.g1_bam,
    g2_bam=rules.SNPsplit_bam.output.g2_bam,
    SNP_bed=rules.get_SNP_bed.output.exons_SNP_bed
  output:
    os.path.join(config["seq_output_dir"], "alignment/STAR/maskedGenome/SNPsplit/bedcov/{sample}_Aligned.out.unique.majchr.dedupl.split.SNP_cov.txt")
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    "{config[script_dir]}sh/allelic_bedcov.sh {threads} {input.unassigned_bam} {input.g1_bam} {input.g2_bam} {input.SNP_bed} {output}"

rule count_genes_from_split_reads:
  input:
    flagged_bam=rules.SNPsplit_bam.output.flagged_bam,
    unassigned_bam=rules.SNPsplit_bam.output.unassigned_bam, 
    g1_bam=rules.SNPsplit_bam.output.g1_bam,
    g2_bam=rules.SNPsplit_bam.output.g2_bam
  output:
    all_counts=os.path.join(config["seq_output_dir"], "featureCounts/SNPsplit/{sample}_all_counts.txt"), 
    unassigned_counts=os.path.join(config["seq_output_dir"], "featureCounts/SNPsplit/{sample}_unassigned_counts.txt"), 
    g1_counts=os.path.join(config["seq_output_dir"], "featureCounts/SNPsplit/{sample}_g1_counts.txt"),
    g2_counts=os.path.join(config["seq_output_dir"], "featureCounts/SNPsplit/{sample}_g2_counts.txt"),
    joined_counts=os.path.join(config["seq_output_dir"], "featureCounts/SNPsplit/{sample}_joined_counts.txt")
  params:
    annotation_file=config["gene_annotation_gtf"],
    gene_counts_param="-F 'GTF' -s 0 -t 'exon' -g 'gene_id' -Q 0 --minOverlap 1 --fracOverlap 0 --fracOverlapFeature 0 -p -C"
  threads: 6
  conda: config["envs_folder"] + "SNPsplit.yml"
  shell:
    """
    featureCounts -a {params.annotation_file} {params.gene_counts_param} -o {output.all_counts} -T {threads} {input.flagged_bam}
    featureCounts -a {params.annotation_file} {params.gene_counts_param} -o {output.unassigned_counts} -T {threads} {input.unassigned_bam}
    featureCounts -a {params.annotation_file} {params.gene_counts_param} -o {output.g1_counts} -T {threads} {input.g1_bam}
    featureCounts -a {params.annotation_file} {params.gene_counts_param} -o {output.g2_counts} -T {threads} {input.g2_bam}
    #merging all featureCounts output in one file
    paste <(awk 'NR>2' {output.all_counts} | sort -k1 | cut -f1,7) <(awk 'NR>2' {output.unassigned_counts} | sort -k1 | cut -f7) <(awk 'NR>2' {output.g1_counts} | sort -k1 | cut -f7) <(awk 'NR>2' {output.g2_counts} | sort -k1 | cut -f7) > {output.joined_counts}
    """

rule verify_trimming:
  input: 
    expand("{path}trimming/{sample}_{mate}.sorted_val_{mate}.fq.gz", path=config["seq_output_dir"], mate=["1","2"], sample=PAIREDSAMPLES)
  output:
    touch(os.path.join(config["logs_dir"], "trimming.done"))

rule verify_STAR:
  input: 
    expand("{path}alignment/STAR/maskedGenome/{sample}_Aligned.out.bam", path=config["seq_output_dir"], sample=PAIREDSAMPLES)
  output:
    touch(os.path.join(config["logs_dir"], "STAR_maskedGenome.done"))

rule verify_SNPsplit:
  input:
    expand("{path}alignment/STAR/maskedGenome/SNPsplit/{sample}_Aligned.out.unique.majchr.nsorted.dedupl.{al_type}.bam", path=config["seq_output_dir"], sample=PAIREDSAMPLES, al_type=["unassigned","genome1","genome2"])
  output:
    touch(os.path.join(config["logs_dir"], "SNPsplit.done"))

rule verify_count_split_reads:
  input:
    expand("{path}featureCounts/SNPsplit/{sample}_joined_counts.txt", path=config["seq_output_dir"], sample=PAIREDSAMPLES)
  output:
    touch(os.path.join(config["logs_dir"], "featureCounts_split_reads.done"))

rule mv_trimmed_fastqc:
  input:
    rules.verify_trimming.output
  output:
    os.path.join(config["logs_dir"], "mv_trimmed_fastqc.done")
  threads: 1
  shell:
    """
    mkdir -p {config[seq_output_dir]}qc/trimmed
    mv {config[seq_output_dir]}trimming/*fastqc* {config[seq_output_dir]}qc/trimmed/ && touch {output}
    """

rule STAR_stats:
  input:
    rules.verify_STAR.output
  output:
    STAR_input=os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/STAR_input_reads.txt"),
    STAR_unique=os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/STAR_uniquely_mapped.txt"),
    STAR_multi=os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/STAR_multiple_loci.txt"),
    STAR_too_many=os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/STAR_too_many_loci.txt")
  threads: 1
  shell:
    """
    grep 'Number of input reads' {config[seq_output_dir]}alignment/STAR/maskedGenome/*_Log.final.out > {output.STAR_input}
    grep 'Uniquely mapped reads number' {config[seq_output_dir]}alignment/STAR/maskedGenome/*_Log.final.out > {output.STAR_unique}
    grep 'Number of reads mapped to multiple loci' {config[seq_output_dir]}alignment/STAR/maskedGenome/*_Log.final.out > {output.STAR_multi}
    grep 'Number of reads mapped to too many loci' {config[seq_output_dir]}alignment/STAR/maskedGenome/*_Log.final.out > {output.STAR_too_many}
    """

rule SNPsplit_stats:
  input:
    rules.verify_SNPsplit.output
  output:
    tot_input_al=os.path.join(config["seq_output_dir"], "stats/SNPsplit/tot_input_al.txt"),
    unassigned=os.path.join(config["seq_output_dir"], "stats/SNPsplit/unassigned_reads.txt"),
    g1_reads=os.path.join(config["seq_output_dir"], "stats/SNPsplit/g1_reads.txt"),
    g2_reads=os.path.join(config["seq_output_dir"], "stats/SNPsplit/g2_reads.txt")
  threads: 1
  shell:
    """
    grep 'read alignments in total' {config[seq_output_dir]}alignment/STAR/maskedGenome/SNPsplit/*SNPsplit_report.txt > {output.tot_input_al}
    grep 'reads were unassignable' {config[seq_output_dir]}alignment/STAR/maskedGenome/SNPsplit/*SNPsplit_report.txt > {output.unassigned}
    grep 'reads were specific for genome 1' {config[seq_output_dir]}alignment/STAR/maskedGenome/SNPsplit/*SNPsplit_report.txt > {output.g1_reads}
    grep 'reads were specific for genome 2' {config[seq_output_dir]}alignment/STAR/maskedGenome/SNPsplit/*SNPsplit_report.txt > {output.g2_reads}
    """

rule spit_reads_counts_stats:
  input:
    rules.verify_count_split_reads.output
  output:
    assigned_all=os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/featureCounts_assigned_all.txt"),
    assigned_unassigned=os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/featureCounts_assigned_un.txt"),
    assigned_g1=os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/featureCounts_assigned_g1.txt"),
    assigned_g2=os.path.join(config["seq_output_dir"], "stats/STAR_SNPsplit/featureCounts_assigned_g2.txt")
  threads: 1
  shell:
    """
    grep 'Assigned' {config[seq_output_dir]}featureCounts/SNPsplit/*_all_counts.txt.summary > {output.assigned_all}
    grep 'Assigned' {config[seq_output_dir]}featureCounts/SNPsplit/*_unassigned_counts.txt.summary > {output.assigned_unassigned}
    grep 'Assigned' {config[seq_output_dir]}featureCounts/SNPsplit/*_g1_counts.txt.summary > {output.assigned_g1}
    grep 'Assigned' {config[seq_output_dir]}featureCounts/SNPsplit/*_g2_counts.txt.summary > {output.assigned_g2}
    """

onsuccess:
  if os.path.isdir('./tmp'):
    shell("rm -r tmp/")
    print("Workflow finished, no error")
  else:
    print("Workflow finished, no error")
