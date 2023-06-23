#!/bin/bash

# RNAseq_regular.sh
#
# A pipeline script for processing RNA-seq data using SLURM workload manager and apptainer container system. 
# The script identifies the fastq files for each sample and initiates a SLURM array to execute the pipeline steps for each sample in parallel. 
# The steps include quality control, trimming, strandedness check, alignment, and quantification.
# See the config_file.sh for module dependencies, version numbers, and reference files used.
#
# List of output files
# - for each sample, trimmed and filtered fastq files (<sample_id>_trimmed_R1_001.fastq.gz and <sample_id>_trimmed_R2_001.fastq.gz)
# - trim_stats.txt: contains number and percentage of total reads after trimming
# - .txt files for each sample indicating read strandedness, as well as a summary .txt file 
# - BAM files containing reads aligned with STAR, as well as BAM index files generated with samtools
# - samtools stats files containing read alignment stats (e.g. unmapped, multi-mapped, etc.)
# - Read counts table produced by trimgalore
# 
# Authors:
# Damon Hofman (d.a.hofman-3@prinsesmaximacentrum.nl)
# Jip van Dinter (j.t.vandinter-3@prinsesmaximacentrum.nl)
#
# Date: 21-06-2023

set -uo pipefail

function usage() {
    cat <<EOF
SYNOPSIS
  run_id_rnaseq_alignment.sh [-c <config file>] [-h]
DESCRIPTION
  run_id the RNAseq processing pipeline consisting of the following steps:
  1. Trim reads with trimgalore (wrapper for cutadapt and fastqc)
  2. Check strandedness of reads with 'howarewestrandedhere'
  3. Map reads with STAR
  4. Create BAM index file and generate read stats with samtools
  5. Quantify gene-level CDS read counts with featureCounts
OPTIONS
  -c, --config <file>    Configuration file to use
  -h, --help             Display this help message
AUTHOR
  Jip van Dinter, MSc
  Damon Hofman, MSc
EOF
}

# Parse command-line arguments
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -c|--config)
    CONFIG="$2"
    shift
    shift
    ;;
    -h|--help)
    usage
    exit
    ;;
    "")
    echo "Error: no option provided"
    usage
    exit 1
    ;;
    *)
    echo "Unknown option: $key"
    usage
    exit 1
    ;;
esac
done

# Check that configuration file is provided
if [[ -z ${CONFIG+x} ]]; then 
    echo "Error: no configuration file provided"
    usage
    exit 1
fi

# Load configuration variables
source $CONFIG

# Load general functions
source ${scriptdir}/general_functions.sh

# Create a unique prefix for the names for this run_id of the pipeline
# This makes sure that run_ids can be identified
run_id=$(uuidgen | tr '-' ' ' | awk '{print $1}')

# Find samples
echo "$(date '+%Y-%m-%d %H:%M:%S') Finding samples..."
get_samples $project_data_folder $data_folder $paired_end

printf "%s\n" "${r1_files[@]}" > ${project_folder}/documentation/r1_files.txt
printf "%s\n" "${r2_files[@]}" > ${project_folder}/documentation/r2_files.txt
printf "%s\n" "${sample_ids[@]}" > ${project_folder}/documentation/sample_ids.txt

# Create output directories
mkdir -p ${project_folder}/log/${run_id}/{trimgalore,star,samtools,howarewestrandedhere} 
mkdir -p ${outdir}

##############################################################################

# Step 1: quality control
trim_jobid=()

trim_jobid+=($(sbatch --parsable \
  --mem=24G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.trimgalore \
  --output=${project_folder}/log/${run_id}/trimgalore/%A_%a.out \
  --export=ALL \
  "${scriptdir}/trimgalore.sh"
))
if [[ ${#trim_jobid[@]} -eq 0 ]]; then
  fatal "TrimGalore job not submitted successfully, trim_jobid array is empty"
fi
info "TrimGalore jobid: ${trim_jobid[@]}"

# Step 2: Check strandedness
strand_jobid=()

strand_jobid+=($(sbatch --parsable \
  --mem=10G \
  --cpus-per-task=2 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.howarewestrandedhere \
  --output=${project_folder}/log/${run_id}/howarewestrandedhere/%A_%a.out \
  --export=ALL \
  "${scriptdir}/check_strand.sh"
))
info "Howarewestrandedhere jobid: ${strand_jobid[@]}"

# Step 3: Align with STAR
star_jobid=()
star_jobid+=($(sbatch --parsable \
  --mem=100G \
  --cpus-per-task=16 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.star_align \
  --output=${project_folder}/log/${run_id}/star/%A_%a.out \
  --dependency=aftercorr:${trim_jobid} \
  --export=ALL \
  "${scriptdir}/star_align.sh"
))
info "STAR alignment jobid: ${star_jobid[@]}"

# Step 4: Create bam index file and generate mapping statistics with samtools
samtools_jobid=()
samtools_jobid+=($(sbatch --parsable \
  --mem=48G \
  --cpus-per-task=16 \
  --gres=tmpspace:100G \
  --time=24:00:00 \
  --job-name=${run_id}.samtools \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --output=${project_folder}/log/${run_id}/samtools/%A_%a.out \
  --dependency=aftercorr:${star_jobid} \
  --export=ALL \
  "${scriptdir}/samtools.sh"
))
info "Samtools jobid: ${samtools_jobid[@]}"
