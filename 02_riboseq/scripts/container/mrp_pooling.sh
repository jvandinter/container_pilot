#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 04-06-2021
#
######################################################################

# Load parameters from main script
threads=$((SLURM_CPUS_PER_TASK * 2))

# Load software modules
module load samtools/${samtools_version}

# Check whether script needs to run
if [[ -f "${outdir}/${pool_id}/${pool_id}.bam" ]]; then
  echo "`date` ${pool_id}.bam already present"
  exit 0
fi

# Create output dirs
mkdir -p "data/processed/${pool_id}"

# Find bam files to merge
bams=$(find ${outdir}/star -maxdepth 2 -name "*.Aligned.sortedByCoord.out.bam" -print)

# Merge bams
echo -e "\n `date` Merging bam files ..."
samtools merge -@ ${threads} "${outdir}/${pool_id}/${pool_id}.bam" $bams
echo -e "\n `date` Merging bam files ... complete! "
