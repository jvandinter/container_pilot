#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 03-06-2021
#
######################################################################

# Load files
mapfile -t r1_files < ${project_folder}/documentation/r1_files.txt
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Set names
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
r1_file="${r1_files[$((SLURM_ARRAY_TASK_ID-1))]}"
r1_filename=$(basename ${r1_file})

# Create output dirs
cd "${outdir}"
mkdir -p "trimgalore/${sample_id}/"

# Check whether script needs to run
if [[ -f "${outdir}/trimgalore/${sample_id}/${r1_filename}" ]]; then
  echo "`date` ${sample_id} file already present"
  exit 0
fi

echo "`date` running trimgalore for ${sample_id}"
apptainer exec -B "/hpc:/hpc" --env "LC_ALL=C.UTF-8" ${container_dir}/trimgalore-0.6.6.sif cutadapt --version
apptainer exec -B "/hpc:/hpc" --env "LC_ALL=C.UTF-8" ${container_dir}/trimgalore-0.6.6.sif fastqc --version
apptainer exec -B "/hpc:/hpc" --env "LC_ALL=C.UTF-8" ${container_dir}/trimgalore-0.6.6.sif trim_galore --version

cd "${outdir}/trimgalore/${sample_id}/"

apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/trimgalore-0.6.6.sif trim_galore \
  "${r1_file}" \
  --cores 2 \
  --gzip \
  --length 25 \
  --trim-n \
  --fastqc \
  --fastqc_args "--outdir ${outdir}/trimgalore/${sample_id}/" \
  --output_dir "${outdir}/trimgalore/${sample_id}/"

mv "${outdir}/trimgalore/${sample_id}/${r1_filename%.*.*}_trimmed.fq.gz" "${outdir}/trimgalore/${sample_id}/${r1_filename}"

# Calculate trimcounts per paired fastq
tot_reads=$(zcat "${r1_file}" | echo $((`wc -l`/4)))
trimmed_reads=$(zcat "${outdir}/trimgalore/${sample_id}/${r1_filename}" | echo $((`wc -l`/4)))
trimmed_percentage=`awk -vn=248 "BEGIN{print(${trimmed_reads}/${tot_reads}*100)}"`

# Add read trimming info to run QC file
printf '%s\t%s\t%s\t%s\n' "${sample_id}_1" "Trimmed" $trimmed_reads $trimmed_percentage >> "${outdir}/trim_stats.txt"

echo "`date` finished ${sample_id}"
