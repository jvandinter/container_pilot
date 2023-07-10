#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 23-06-2023
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
source ${scriptdir}/mrp_functions.sh

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

mkdir -p ${project_folder}/log/${run_id}/{trimgalore,bowtie2,star_align,riboseqc} data/processed

echo -e "\n ====== `date` Map Riboseq Pipeline ====== \n"

echo -e "\n`date` Filtering and trimming ..."
echo -e "====================================================================================== \n"

# 1. TRIMGALORE. Parallel start of all trimgalore jobs to filter for quality
#                with CUTADAPT and output quality reports with FASTQC

trim_jobid=()

trim_jobid+=($(sbatch --parsable \
  --mem=8G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.trimgalore \
  --output=${project_folder}/log/${run_id}/trimgalore/%A_%a.out \
  --export=ALL \
  ${scriptdir}/mrp_trimgalore.sh 
))

if [[ ${#trim_jobid[@]} -eq 0 ]]; then
  fatal "TrimGalore job not submitted successfully, trim_jobid array is empty"
fi

info "trimgalore jobid: ${trim_jobid}"

echo -e "\n`date` Removing contaminants ..."
echo -e "====================================================================================== \n"

# 2. BOWTIE2. Use combination of tRNA, rRNA, snRNA, snoRNA, mtDNA fastas to
#             remove those contaminants from RIBO-seq data. Outputs QC stats
#             to a file per contaminant group.

contaminant_jobid=()

contaminant_jobid+=($(sbatch --parsable \
  --mem=8G \
  --cpus-per-task=12 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.contaminants \
  --output=${project_folder}/log/${run_id}/bowtie2/%A_%a.out \
  --dependency=aftercorr:${trim_jobid} \
  --export=ALL \
  ${scriptdir}/mrp_remove_contaminants.sh
  
))

info "contaminant jobid: ${contaminant_jobid}"

echo -e "\n`date` Align reads to genome with STAR ..."
echo -e "====================================================================================== \n"

# 3. STAR. Align contaminant-depleted read files to supplied genome and
#          transcriptome. If no new custom transcriptome is supplied, 
#          the normal reference transcriptome is used for 
#          guided assembly.

star_jobid=()

star_jobid+=($(sbatch --parsable \
  --mem=80G \
  --cpus-per-task=12 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.star_align \
  --output=${project_folder}/log/${run_id}/star_align/%A_%a.out \
  --dependency=aftercorr:${contaminant_jobid} \
  --export=ALL \
  ${scriptdir}/mrp_star_align.sh
))

info "alignment jobid: ${star_jobid}"

echo -e "\n`date` Perform QC with RiboseQC ..."
echo -e "====================================================================================== \n"

# 4. RiboseQC. 

riboseqc_jobid=()

riboseqc_jobid+=($(sbatch --parsable \
  --mem=4G \
  --cpus-per-task=1 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.riboseqc \
  --output=${project_folder}/log/${run_id}/riboseqc/%A_%a.out \
  --dependency=aftercorr:${star_jobid} \
  --export=ALL \
  ${scriptdir}/mrp_riboseqc.sh
))

info "RiboseQC jobid: ${riboseqc_jobid}"

echo -e "\n`date` Creating MultiQC reports ..."
echo -e "====================================================================================== \n"

# 5. MultiQC.

# multiqc_jobid=()

# multiqc_jobid+=($(sbatch --parsable \
#   --mem=${low_mem} \
#   --cpus-per-task=${low_cpu} \
#   --time=24:00:00 \
#   --job-name=${run_id}.multiqc \
#   --output=${project_folder}/log/${run_id}/%A_multiqc.out \
#   --dependency=afterok:${riboseqc_jobid} \
#   --export=ALL \
#   ${scriptdir}/mrp_multiqc.sh
# ))

# info "MultiQC jobid: ${multiqc_jobid[@]}"

echo -e "\n ====== `date` Started all jobs! ====== \n"
