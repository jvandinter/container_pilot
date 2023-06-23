#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 03-06-2021
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
if [[ -f "${wd}/data/processed/RiboseQC/${sample_id}/${sample_id}.html" ]]; then
  echo "`date` ${sample_id} QC report already present"
  exit 0
fi

echo "`date` doing RIBO-seq QC for ${sample_id}"

# Create output dirs
cd "${wd}/data/processed"
mkdir -p "RiboseQC/${sample_id}/"

# Use RiboSeQC to generate HTML report of the data
Rscript "${scriptdir}/mrp_riboseqc.R" \
  "${wd}/data/processed/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam" \
  "${wd}/data/processed/RiboseQC/${sample_id}/${sample_id}" \
  "${rannot}" \
  "${annot_name}" \
  "${pandoc_dir}" \
  "${resource_dir}" \
  "${annotation_package}"

echo "`date` finished ${sample_id}"
