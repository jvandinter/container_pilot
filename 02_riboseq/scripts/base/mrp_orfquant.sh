#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 04-06-2021
#
######################################################################

# Load parameters from main script
source $1
rannot=$2
annotation_package=$3
cpu=$4
threads=$((cpu * 2))

# Get sample IDs
fastq_files=(${wd}/data/raw/*fastq.gz)

# Initiate arrays
    sample_ids=()

# Get sample IDs and barcodes from fastq files
for i in ${!fastq_files[@]}; do

    sample_ids[i]=$(basename ${fastq_files[i]} | cut -f 1 -d "_")

done

sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Load software modules
module load R/${r_version}

# Check whether script needs to run
if [[ -f "${wd}/data/processed/ORFquant/${sample_id}/${sample_id}.html" ]]; then
  echo "File already present"
  exit 0
fi

echo "`date` calling ORFs for ${sample_id}"

# Create output dirs
cd "${wd}/data/processed/"
mkdir -p "ORFquant/${sample_id}/"

Rscript "${scriptdir}/mrp_orfquant.R" \
  ${wd} \
  "${wd}/data/processed/RiboseQC/${sample_id}/${sample_id}_for_ORFquant" \
  "${sample_id}" \
  "${annot_name}" \
  "${rannot}" \
  "${threads}" \
  "${pandoc_dir}" \
  "${resource_dir}" \
  "${annotation_package}"

  echo "`date` finished ${sample_id}"
