#!/bin/bash

# run_rnaseq_transcriptome_assembly.sh
#
# Short description of pipeline script
#
#

set -uo pipefail

function usage() {
    cat <<EOF
SYNOPSIS
  run_rnaseq_transcriptome_assembly.sh [-c <config file>] [-h]
DESCRIPTION
  1. Assemble transcriptomes for each sample with STRINGTIE, using BAM files as input
  2. Run GFFCOMPARE to compare custom GTF to reference GTF 
  3. Merge GTF files per sample group with STRINGTIE --merge
  4. Annotate and filter merged GTF with GFFCOMPARE and custom R scripts
  5. Create custom annotation with custom merged GTF for ribo-seq analysis
  6. Generate Salmon index
  7. Quantify reads with salmon quant
  8. Run MultiQC on new data
OPTIONS
  -c, --config <file>    Configuration file to use
  -h, --help             Display this help message
AUTHOR
  Jip van Dinter, MSc
  Damon Hofman, MSc
EOF
}

# Parse command-line arguments
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -c|--config)
    CONFIG="$2"
    shift
    shift
    ;;
    -h|--help)
    usage
    exit
    ;;
    "")
    echo "Error: no option provided"
    usage
    exit 1
    ;;
    *)
    echo "Unknown option: $key"
    usage
    exit 1
    ;;
esac
done

# Check that configuration file is provided
if [[ -z ${CONFIG+x} ]]; then 
    echo "Error: no configuration file provided"
    usage
    exit 1
fi

# Load configuration variables
source $CONFIG

# Load general functions
source ${scriptdir}/general_functions.sh

# Find samples
echo "$(date '+%Y-%m-%d %H:%M:%S') Finding samples..."
get_samples $project_data_folder $data_folder $paired_end

# Create a unique prefix for the names for this run of the pipeline. 
# This makes sure that runs can be identified
run=$(uuidgen | tr '-' ' ' | awk '{print $1}')

# Create output directories
mkdir -p log/${run}/{gffcompare,stringtie,salmon_quant}

################################################################################

# Step 1: Transcript assembly with stringtie
stringtie_jobid=()
stringtie_jobid+=($(sbatch --parsable \
  --mem=24G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.stringtie \
  --output=log/${run}/stringtie/%A_%a.out \
  --export=ALL\
  ${scriptdir}/stringtie.sh
))
info "stringtie jobid: ${stringtie_jobid[@]}"

# Step 2: Compare sample GTF with reference GTF using gffcompare
gffcompare_jobid=()
gffcompare_jobid+=($(sbatch --parsable \
  --mem=10G \
  --cpus-per-task=2 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.gffcompare \
  --output=log/${run}/gffcompare/%A_%a \
  --dependency=aftercorr:${stringtie_jobid} \
  --export=ALL \
  ${scriptdir}/gffcompare.sh
))
info "Gffcompare jobid: ${gffcompare_jobid[@]}"

# Step 3: Merge all sample GTF files into single non-redundant GTF file using reference GTF file as guide
stringtie_merge_jobid=()
stringtie_merge_jobid+=($(sbatch --parsable \
  --mem=24G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --job-name=${run}.stringtie_merge \
  --output=log/${run}/%A_stringtie_merge.out \
  --dependency=afterany:${stringtie_jobid} \
  --export=ALL \
  ${scriptdir}/stringtie_merge.sh
))
info "Stringtie merge jobid: ${stringtie_merge_jobid[@]}"

# Step 4: Annotate transcripts with gffcompare class codes and perform transcript filtering
filter_annotate_jobid=()
filter_annotate_jobid+=($(sbatch --parsable \
  --mem=24G \
  --cpus-per-task=2 \
  --time=24:00:00 \
  --job-name=${run}.filter_annotate \
  --output=log/${run}/%A_filter_annotate.out \
  --dependency=afterany:${stringtie_merge_jobid} \
  --export=ALL \
  ${scriptdir}/filter_annotate.sh
))
info "Filter and annotation jobid: ${filter_annotate_jobid[@]}"

# Step 5: Create custom annotation file for RiboseQC and ORFquant based on filtered custom GTF
custom_annotation_jobid=()
if [[ ${create_annotation} =~ "TRUE" ]]; then
  if [[ $(find ${wd}/ -name '*.Rannot' | wc -l) -eq 0 ]]; then
    custom_annotation_jobid+=($(sbatch --parsable \
      --mem=10G \
      --cpus-per-task=2 \
      --gres=tmpspace:50G \
      --time=24:00:00 \
      --job-name=${run}.custom_annotation \
      --output=log/${run}/%A_custom_annotation.out \
      --dependency=afterok:${filter_annotate_jobid} \
      --export=ALL \
      ${scriptdir}/custom_annotation.sh
    ))
    info "Custom annotation jobid: ${custom_annotation_jobid[@]}"
  else
  echo "Annotation file already present"
  fi
else
  echo "Creation of annotation not specified"
fi