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
module load R/${r_version}

# Load files
mapfile -t r1_files < ${project_folder}/documentation/r1_files.txt
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Set names
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Check whether script needs to run
if [[ -f "${outdir}/ORFquant/${sample_id}/${sample_id}.html" ]]; then
  echo "File already present"
  exit 0
fi

echo "`date` calling ORFs for ${sample_id}"

# Create output dirs
cd "${outdir}/"
mkdir -p "ORFquant/${sample_id}/"

apptainer exec -B "/hpc:/hpc" ${container_dir}/orfquant-4.1.2.sif \ 
  Rscript "${scriptdir}/mrp_orfquant.R" \
  ${wd} \
  "${outdir}/RiboseQC/${sample_id}/${sample_id}_for_ORFquant" \
  "${sample_id}" \
  "${annot_name}" \
  "${rannot}" \
  "${threads}" \
  "${pandoc_dir}" \
  "${resource_dir}" \
  "${annotation_package}"  \
  ${package_install_loc}

  echo "`date` finished ${sample_id}"
