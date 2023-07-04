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
mapfile -t r1_files < ${project_folder}/documentation/r1_files.txt
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Set names
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
r1_file="${r1_files[$((SLURM_ARRAY_TASK_ID-1))]}"
r1_filename=$(basename ${r1_file})

# Check whether script needs to run
if [[ -s "${outdir}/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam.bai" ]]; then
  echo "`date` ${sample_id} already present"
  exit 0
fi

echo "`date` mapping ${sample_id}"

# Create output dirs
cd "${outdir}"
mkdir -p "star/${sample_id}/"

# Map ribo reads
apptainer exec -B /hpc:/hpc ${container_dir}/STAR-2.7.8a.sif STAR --version
apptainer exec -B /hpc:/hpc ${container_dir}/STAR-2.7.8a.sif STAR \
  --genomeDir "${star_index_basedir}/29nt" \
  --sjdbGTFfile ${gtf} \
  --runThreadN ${threads} \
  --runDirPerm All_RWX \
  --twopassMode Basic \
  --readFilesIn \
  "${outdir}/bowtie2/${sample_id}/${sample_id}_filtered.fastq.gz" \
  --readFilesCommand zcat \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 20 \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix "${outdir}/star/${sample_id}/${sample_id}." \
  --limitOutSJcollapsed 10000000 \
  --limitIObufferSize=300000000 \
  --outFilterType BySJout \
  --alignSJoverhangMin 1000 \
  --outTmpKeep None

echo "`date` indexing ${sample_id}"

apptainer exec -B /hpc:/hpc ${container_dir}/samtools-1.12.sif samtools index \
  -@ ${threads} \
  "${outdir}/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam"

echo "`date` finished ${sample_id}"
