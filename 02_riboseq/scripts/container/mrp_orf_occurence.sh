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

# Load correct modules
module load R/${r_version}

# Check whether script needs to run
if [[ -f "${wd}/analysis/${pool_id}/${pool_id}_ORFcounts.csv" ]]; then
  echo "File already present"
  exit 0
fi

orfquant_results="${wd}/analysis/${pool_id}/ORFquant/${pool_id}_final_ORFquant_results"

# Count ORFs across samples and export table with pooled stats
echo -e "\n `date` Counting ORFs across samples ..."
Rscript "${scriptdir}/mrp_orf_occurence.R" \
  ${wd} \
  ${orfquant_results} \
  ${pool_id}