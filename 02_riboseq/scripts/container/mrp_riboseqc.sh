#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 03-06-2021
#
######################################################################

# Load parameters from main script

threads=$((SLURM_CPUS_PER_TASK * 2))

# Load files
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Set names
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Check whether script needs to run
if [[ -f "${outdir}/RiboseQC/${sample_id}/${sample_id}.html" ]]; then
  echo "`date` ${sample_id} QC report already present"
  exit 0
fi

echo "`date` doing RIBO-seq QC for ${sample_id}"

# Create output dirs
cd "${outdir}"
mkdir -p "RiboseQC/${sample_id}/"

# Use RiboSeQC to generate HTML report of the data
apptainer exec -B "/hpc:/hpc" ${container_dir}/orfquant-4.1.2.sif \
  Rscript "${scriptdir}/mrp_riboseqc.R" \
  "${outdir}/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam" \
  "${outdir}/RiboseQC/${sample_id}/${sample_id}" \
  "${rannot}" \
  "${annot_name}" \
  "${pandoc_dir}" \
  "${resource_dir}" \
  "${annotation_package}"

echo "`date` finished ${sample_id}"
