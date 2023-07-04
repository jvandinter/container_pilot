#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 04-06-2021
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
if [[ -f "${wd}/data/processed/ORFquant/${pool_id}/${pool_id}.html" ]]; then
  echo "File already present"
  exit 0
fi

echo "`date` calling ORFs for ${pool_id}"

# Create output dirs
cd "${wd}/data/processed/"
mkdir -p "ORFquant/${pool_id}/"

Rscript "${scriptdir}/mrp_orfquant.R" \
  ${wd} \
  "${wd}/data/processed/RiboseQC/${pool_id}/${pool_id}_for_ORFquant" \
  "${pool_id}" \
  "${annot_name}" \
  "${rannot}" \
  "${threads}" \
  "${pandoc_dir}" \
  "${resource_dir}" \
  "${annotation_package}"

  echo "`date` finished ${pool_id}"
