#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 03-06-2021
#
######################################################################

# Load correct modules
module load cutadapt/${cutadapt_version}
module load fastqc/${fastqc_version}
module load trimgalore/${trimgalore_version}

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
if [[ -f "${outdir}/trimgalore/${sample_id}/${r1_filename/_R1_/_R1_trimmed_}" ]]; then
  echo "`date` ${sample_id} file already present"
  exit 0
fi

echo "`date` running trimgalore for ${sample_id}"
cutadapt --version
fastqc --version

cd "${outdir}/trimgalore/${sample_id}/"

trim_galore "${r1_file}" \
  --cores 2 \
  --gzip \
  --length 25 \
  --trim-n \
  --fastqc \
  --fastqc_args "--outdir ${wd}/data/processed/trimgalore/${sample_id}/" \
  --output_dir "$wd//data/processed/trimgalore/${sample_id}/"

mv "${outdir}/trimgalore/${sample_id}/${fastq%.*.*}_trimmed.fq.gz" "${outdir}/trimgalore/${sample_id}/${fastq}"

# Calculate trimcounts per paired fastq
tot_reads=$(zcat "${r1_file}" | echo $((`wc -l`/4)))
trimmed_reads=$(zcat "${outdir}/trimgalore/${sample_id}/${r1_trimmed}" | echo $((`wc -l`/4)))
trimmed_percentage=`awk -vn=248 "BEGIN{print(${trimmed_reads}/${tot_reads}*100)}"`

# Add read trimming info to run QC file
printf '%s\t%s\t%s\t%s\n' "${sample_id}_1" "Trimmed" $trimmed_reads $trimmed_percentage >> "${outdir}/trim_stats.txt"

echo "`date` finished ${sample_id}"
