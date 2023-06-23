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
gtf=$2
cpu=$3
threads=$((cpu * 2))

# Load correct modules
module load STAR/${star_version}
module load samtools/${samtools_version}

# Get sample IDs
fastq_files=(${wd}/data/raw/*fastq.gz)

# Initiate arrays
    sample_ids=()

# Get sample IDs and barcodes from fastq files
for i in ${!fastq_files[@]}; do

    sample_ids[i]=$(basename ${fastq_files[i]} | cut -f 1 -d "_")

done

sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Check whether script needs to run
if [[ -s "${wd}/data/processed/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam.bai" ]]; then
  echo "`date` ${sample_id} already present"
  exit 0
fi

echo "`date` mapping ${sample_id}"

# Create output dirs
cd "${wd}/processed"
mkdir -p "star/${sample_id}/"

# Map ribo reads
STAR --genomeDir "${star_index_basedir}/29nt" \
  --sjdbGTFfile ${gtf} \
  --runThreadN ${threads} \
  --runDirPerm All_RWX \
  --twopassMode Basic \
  --readFilesIn \
  "${wd}/data/processed/bowtie2/${sample_id}/${sample_id}_filtered.fastq.gz" \
  --readFilesCommand zcat \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 20 \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix "${wd}/data/processed/star/${sample_id}/${sample_id}." \
  --limitOutSJcollapsed 10000000 \
  --limitIObufferSize=300000000 \
  --outFilterType BySJout \
  --alignSJoverhangMin 1000 \
  --outTmpKeep None

echo "`date` indexing ${sample_id}"

samtools index -@ ${threads} "${wd}/data/processed/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam"

echo "`date` finished ${sample_id}"
