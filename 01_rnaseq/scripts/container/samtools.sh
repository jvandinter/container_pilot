#!/bin/bash

set -uo pipefail

# Get correct files
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
bam="${sample_id}.Aligned.out.bam"
new_bam="${sample_id}.Aligned.sortedByCoord.out.bam"

echo "`date` running samtools for ${sample_id}"

# Check whether script needs to run
if [[ -f "${outdir}/star/${sample_id}/${new_bam}.bai" ]]; then
  echo "`date` ${sample_id} BAM index already present"
  exit 0
fi

if ! [[ -s "${outdir}/star/${sample_id}/${bam}" ]]; then
  echo "${sample_id} BAM file is empty. Re-run STAR."
  exit 1
fi

# Sort BAM
apptainer exec -B /hpc:/hpc,${TMPDIR}:${TMPDIR} ${container_dir}/samtools-1.12.sif samtools sort \
  -@ 8 \
  -l 9 \
  -o "${outdir}/star/${sample_id}/${new_bam}" \
  -T "${TMPDIR}" \
  "${outdir}/star/${sample_id}/${bam}"

# Create mapping statistics with samtools
apptainer exec -B /hpc:/hpc ${container_dir}/samtools-1.12.sif samtools stats -@ 8 \
  "${outdir}/star/${sample_id}/${new_bam}" > "${outdir}/star/${sample_id}/${sample_id}_stats.txt"

# Index the bam with samtools
apptainer exec -B /hpc:/hpc ${container_dir}/samtools-1.12.sif samtools index -@ 8 "${outdir}/star/${sample_id}/${new_bam}"

# Remove original STAR bam
if [[ -s "${outdir}/star/${sample_id}/${new_bam}" ]]; then
  rm "${outdir}/star/${sample_id}/${bam}"
fi

echo "`date` finished indexing ${sample_id}"
