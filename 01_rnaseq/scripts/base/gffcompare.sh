#!/bin/bash

set -uo pipefail

# Load modules
module load gffcompare/${gffcompare_version}

# Load sample ids
mapfile -t sample_ids < sample_ids.txt
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Check whether script needs to run
if [[ -f "${outdir}/gffcompare/${sample_id}/${sample_id}.stats" ]]; then
  echo "`date` gffcompare already completed for ${sample_id}"
  exit 0
fi

# Create output dirs
mkdir -p "${outdir}/gffcompare/${sample_id}"
cd "${outdir}/gffcompare/${sample_id}/"

# Run GFFcompare to compare sample GTF to reference GTF
gffcompare --version
gffcompare \
  -r ${reference_gtf} \
  -o "${sample_id}" \
  "${outdir}/stringtie/${sample_id}/${sample_id}.gtf"

echo "`date` GFFCompare for ${sample_id} finished!"