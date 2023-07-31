#!/bin/bash

set -uo pipefail

# Check whether script needs to run
if [[ -f "${outdir}/customannotation/${merged_gtf_basename}/${merged_gtf_basename}_novel_filtered.gtf_Rannot" ]]; then
  echo "Custom annotation package already generated!"
  exit 0
fi

echo "`date` Generating custom annotation ..."

# Export environment for apptainer

# Prepare RiboseQC and ORFquant annotation files based on 
# merged annotated stringtie GTF
apptainer exec -B "/hpc:/hpc",${TMPDIR}:${TMPDIR} --env "LC_CTYPE=en_US.UTF-8" ${container_dir}/orfquant-4.1.2.sif \
  Rscript "${scriptdir}/prepare_custom_annotation.R" \
  ${twobit} \
  "${outdir}/customannotation/${merged_gtf_basename}_novel_filtered.gtf" \
  "${outdir}/customannotation/${merged_gtf_basename}/" \
  ${merged_gtf_basename} \
  ${package_install_loc}

echo "`date` Generating custom annotation finished"
