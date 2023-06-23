#!/bin/bash

#SBATCH -t 4:00:00
#SBATCH --job-name=ribo-sebas

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 11-10-2022
#
######################################################################


function usage() {
    cat <<EOF
SYNOPSIS
  mrp_main_script.sh ./path/to/config.config - run orf identifier pipeline
  mrp_main_script.sh help - display this help message
DESCRIPTION
  01. Pool BAMs with SAMTOOLS
  02. Run pooled RIBOSEQC
  03. Run pooled ORFQUANT

AUTHOR
  Jip van Dinter, MSc
  Damon Hofman, MSc
EOF
}

function info() {
    echo "INFO: $@" >&2
}
function error() {
    echo "ERR:  $@" >&2
}
function fatal() {
    echo "ERR:  $@" >&2
    exit 1
}

# Create a unique prefix for the names for this run of the pipeline. 
# This makes sure that runs can be identified
run=$(uuidgen | tr '-' ' ' | awk '{print $1}')

echo "`date` run ID $run"

# Show help message if there was no config file location given on the commandline
if [[ -z $1 ]]; then 

  usage; exit;

fi

# Source all variables from the config file
CONFIG=$1
source ${CONFIG}
source ${scriptdir}/mrp_functions.sh

################################################################################
#
# Find fastq samples in directory
#
################################################################################

# Get sample IDs
fastq_files=(${wd}/data/raw/*fastq.gz)

# Initiate arrays
    sample_ids=()
    samples=()

# Get sample IDs and barcodes from fastq files
for i in ${!fastq_files[@]}; do

    sample_ids[i]=$(basename ${fastq_files[i]} | cut -f 1 -d "_")
    samples[i]=$(basename ${fastq_files[i]})

done

check_annotation ${reference_annotation} ${reference_gtf} ${custom_annotation} ${custom_gtf} ${custom_annotation_package} ${reference_annotation_package}

echo "`date` using ${annotation_package}"
echo "`date` using ${rannot}"

# make sure there are samples
if [[ ${#samples[@]} -eq 0 ]]; then
  fatal "no samples found in ./raw/ or file containing fastq file locations not present"
fi

info "samples: n = ${#samples[@]}"
for sample in ${samples[@]}; do
  info "    $sample"
done

################################################################################
#
# Run the pipeline
#
################################################################################

mkdir -p log/${run}/

echo -e "\n ====== `date` Map Riboseq Pipeline ====== \n"

echo -e "\n`date` Merge BAMs for pooled analysis ..."
echo -e "====================================================================================== \n"

# 1. SAMTOOLS. Pooling of samples 

pool_jobid=()

pool_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${medium_cpu} \
  --time=24:00:00 \
  --job-name=${run}.pooling \
  --output=log/${run}/pooling.out \
  ${scriptdir}/mrp_pooling.sh \
  ${CONFIG} \
  ${medium_cpu}
))

info "Pooling jobid: ${pool_jobid}"

echo -e "\n`date` Perform QC with RiboseQC ..."
echo -e "====================================================================================== \n"

# 2. RiboseQC. Merged RiboseQC

riboseqc_merged_jobid=()

riboseqc_merged_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${medium_cpu} \
  --time=24:00:00 \
  --job-name=${run}.riboseqc.pooled \
  --output=log/${run}/riboseqc.pooled.out \
  --dependency=afterok:${pool_jobid} \
  ${scriptdir}/mrp_merged_riboseqc.sh \
  ${CONFIG} \
  ${rannot} \
  ${annotation_package} \
  ${medium_cpu}
))

info "RiboseQC pooled jobid: ${riboseqc_merged_jobid}"

echo -e "\n`date` Detect ORFs with ORFquant using pooled samples ..."
echo -e "====================================================================================== \n"

# 3. ORFquant. Merged ORFquant

orfquant_merged_jobid=()

orfquant_merged_jobid+=($(sbatch --parsable \
  --mem=${high_mem} \
  --cpus-per-task=${high_cpu} \
  --time=144:00:00 \
  --job-name=${run}.orfquant.pooled \
  --output=log/${run}/orfquant.pooled.out \
  --dependency=afterok:${riboseqc_merged_jobid} \
  ${scriptdir}/mrp_merged_orfquant.sh \
  ${CONFIG} \
  ${rannot} \
  ${annotation_package} \
  ${high_cpu}
))

info "Merged ORFquant jobid: ${orfquant_merged_jobid}"

echo -e "\n ====== `date` Started all jobs! ====== \n"
