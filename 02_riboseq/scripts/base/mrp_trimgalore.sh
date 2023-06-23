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
source $2

# Load correct modules
module load cutadapt/${cutadapt_version}
module load fastqc/${fastqc_version}
module load trimgalore/${trimgalore_version}

# Get correct files
get_samples ${wd} ${data_folder}

full_path_fastq="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
fastq=$(basename ${full_path_fastq})

# Check whether script needs to run
if [[ -f "${wd}/data/processed/trimgalore/${sample_id}/${fastq}" ]]; then
  echo "${sample_id} already trimmed"
  exit 0
fi

echo "`date` running trimgalore for ${sample_id}"

# Create output dirs
cd "${wd}/data/processed"
mkdir -p "trimgalore/${sample_id}/"

# Run trimgalore on both reads
cutadapt --version
fastqc --version

# Change names of trimgalore output
cd "$wd/data/processed/trimgalore/${sample_id}/"

trim_galore "${full_path_fastq}" \
  --cores 2 \
  --gzip \
  --length 25 \
  --trim-n \
  --fastqc \
  --fastqc_args "--outdir ${wd}/data/processed/trimgalore/${sample_id}/" \
  --output_dir "$wd//data/processed/trimgalore/${sample_id}/"

mv "${wd}/data/processed/trimgalore/${sample_id}/${fastq%.*.*}_trimmed.fq.gz" "${wd}/data/processed/trimgalore/${sample_id}/${fastq}"

# Calculate trimcounts per paired fastq
tot_reads_1=$(zcat "${full_path_fastq}" | echo $((`wc -l`/4)))

trimmed_reads_1=$(zcat "${wd}/data/processed/trimgalore/${sample_id}/${fastq}" | echo $((`wc -l`/4)))
trimmed_percentage_1=`awk -vn=248 "BEGIN{print(${trimmed_reads_1}/${tot_reads_1}*100)}"`

# Add read trimming info to run QC file
printf '%s\t%s\t%s\t%s\n' "${sample_id}_1" "Trimmed" $trimmed_reads_1 $trimmed_percentage_1 >> "$wd/data/processed/trim_stats.txt"

echo "`date` finished ${sample_id}"
