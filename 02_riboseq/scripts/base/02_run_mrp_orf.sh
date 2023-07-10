#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 04-07-2023
#
######################################################################

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
source ${scriptdir}//mrp_functions.sh

# Create a unique prefix for the names for this run_id of the pipeline
# This makes sure that run_ids can be identified
run_id=$(uuidgen | tr '-' ' ' | awk '{print $1}')

################################################################################
#
# Find fastq samples in directory
#
################################################################################

# Find samples
echo "$(date '+%Y-%m-%d %H:%M:%S') Finding samples..."
get_samples $project_data_folder $data_folder

printf "%s\n" "${r1_files[@]}" > ${project_folder}/documentation/r1_files.txt
printf "%s\n" "${sample_ids[@]}" > ${project_folder}/documentation/sample_ids.txt

echo "$(date '+%Y-%m-%d %H:%M:%S') Finding Annotation..."
check_annotation ${reference_annotation} ${reference_gtf} ${reference_annotation_package} ${custom_annotation} ${custom_gtf} ${custom_annotation_package}

echo "`date` using ${annotation_package}"
echo "`date` using ${rannot}"
echo "`date` using ${gtf}"

export annotation_package=${annotation_package}
export rannot=${rannot}
export gtf=${gtf}

# Create output directories
mkdir -p ${project_folder}/log/${run_id}/{trimgalore,star_align,bowtie2,riboseqc} 
echo "`date` using run ID: ${run_id}"
mkdir -p ${outdir}

# make sure there are samples
if [[ ${#samples[@]} -eq 0 ]]; then
  fatal "no samples found in ./raw/ or file containing fastq file locations not present"
fi

info "samples: n = ${#samples[@]}"
for sample in ${samples[@]}; do
  info "    $sample"
done

################################################################################
#
# Run the pipeline
#
################################################################################

mkdir -p ${project_folder}/log/${run_id}/orfquant/

echo -e "\n ====== `date` Map Riboseq Pipeline ====== \n"

echo -e "\n`date` Run ORFquant ..."
echo -e "====================================================================================== \n"

orfquant_jobid=()

orfquant_jobid+=($(sbatch --parsable \
  --mem=${high_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.orfquant \
  --output=${project_folder}/log/${run_id}/orfquant/%A_%a.out \
  --export=ALL \
  ${scriptdir}/mrp_orfquant.sh
))

if [[ ${#orfquant_jobid[@]} -eq 0 ]]; then
  fatal "ORFquant job not submitted successfully, trim_jobid array is empty"
fi

info "orfquant jobid: ${orfquant_jobid}"
