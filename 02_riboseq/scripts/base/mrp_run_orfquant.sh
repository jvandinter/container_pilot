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
  02. Run pooled orfquant
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

mkdir -p log/${run}/orfquant/

echo -e "\n ====== `date` Map Riboseq Pipeline ====== \n"

echo -e "\n`date` Run ORFquant ..."
echo -e "====================================================================================== \n"

orfquant_jobid=()

orfquant_jobid+=($(sbatch --parsable \
  --mem=${high_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.orfquant \
  --output=log/${run}/orfquant/%A_%a.out \
  ${scriptdir}/mrp_orfquant.sh \
  ${CONFIG} \
  ${rannot} \
  ${annotation_package} \
  ${high_cpu}
))

info "orfquant jobid: ${orfquant_jobid}"
