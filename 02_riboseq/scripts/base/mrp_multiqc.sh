#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

set -uo pipefail

# Load parameters from main script
source $1
name=$2

# Load correct modules
module load multiqc/${multiqc_version}

echo "`date` running MultiQC for all samples"

cd ${wd}

# Run MultiQC
multiqc \
  "${wd}/data/processed/" \
  --outdir "${wd}/data/processed/" \
  --filename "${pool_id}_multiqc.html"

  echo "`date` finished MultiQC"
