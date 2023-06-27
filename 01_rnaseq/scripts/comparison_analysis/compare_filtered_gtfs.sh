#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem=24G 
#SBATCH --time=24:00:00 

set -uo pipefail

source "/hpc/pmc_vanheesch/projects/Jip/pilots/20230606_JD_containers/01_rnaseq/documentation/container_pilot.config"

container_gtf="${project_folder}/analysis/container/customannotation/container_pancreas_novel_filtered.gtf"
base_gtf="${project_folder}/analysis/base/customannotation/base_merged_novel_filtered.gtf"

mkdir -p "${project_folder}/analysis/gffcompare/pipeline_comparison/"

# Create output dirs
cd "${project_folder}/analysis/gffcompare/pipeline_comparison/"

# Run GFFcompare to annotate novel assembly .GTF
apptainer exec -B /hpc:/hpc ${container_dir}/gffcompare-0.12.6.sif gffcompare \
  -V \
  -r "${base_gtf}" \
  -s "${masked_fasta}" \
  -o "base_vs_container_comparison_gffcompare" \
  "${container_gtf}"

  # Run GFFcompare to annotate novel assembly .GTF
apptainer exec -B /hpc:/hpc ${container_dir}/gffcompare-0.12.6.sif gffcompare \
  -V \
  -r "${container_gtf}" \
  -s "${masked_fasta}" \
  -o "container_vs_base_comparison_gffcompare" \
  "${base_gtf}"

echo "`date` running GFFcompare for transcript occurence"
