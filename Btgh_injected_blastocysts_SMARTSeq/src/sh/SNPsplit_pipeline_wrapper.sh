############################
# Wrapper to run Snakemake #
############################

# The following options are included in every command because they are generally useful
generalOptions="--rerun-incomplete --latency-wait 120 --keep-going"

## explanation: 
## Re-run all jobs the output of which is recognized as incomplete, wait this time before considering a job does not have output, keep going with other jobs if one fails

# Options helpful for debugging
debuggingOptions="--printshellcmds --verbose --reason"

## explanation:
## print out the shell commands that will be executed, the debugging output and the reason for each rule to be executed

# Going to the main folder of the experiment folder, which will be the parent of all folders in the pipeline and in the Snakefile
cd /g/boulard/Projects/Btgh_injected_blastocysts_SMARTSeq/

snakemake \
  -s snake-make/SNPsplit.Snakefile \
  --configfile config/SNPsplit_config.yaml \
  --use-conda \
  $generalOptions\
  $debuggingOptions\
  -j 248 \
  --cluster-config config/SNPsplit_cluster.SLURM.json \
  --cluster "sbatch -J {cluster.name} --cpus-per-task {cluster.nCPUs} --mem {cluster.memory} \
                   --time {cluster.maxTime} -o {cluster.output} -e {cluster.error}" \
  --local-cores 1 \
  --max-jobs-per-second 1

