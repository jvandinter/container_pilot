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
rannot=$2
annotation_package=$3
cpu=$4
threads=$((cpu * 2))

# Load software modules
module load R/${r_version}

# Check whether script needs to run
if [[ -f "${wd}/data/processed/RiboseQC/${pool_id}/${pool_id}.html" ]]; then
  echo "File already present"
  exit 0
fi

echo "`date` generating P sites for ${pool_id}"

# Create output dirs
cd "${wd}/data/processed/" 
mkdir -p "RiboseQC/${pool_id}/"

# Use RiboSeQC to generate HTML report of the data
Rscript "${scriptdir}/mrp_riboseqc.R" \
  "${wd}/data/processed/${pool_id}/${pool_id}.bam" \
  "${wd}/data/processed/RiboseQC/${pool_id}/${pool_id}" \
  "${rannot}" \
  "${annot_name}" \
  "${pandoc_dir}" \
  "${resource_dir}" \
  "${annotation_package}"

echo "`date` finished ${pool_id}"
