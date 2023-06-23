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
source $2
cpu=$3
threads=$((cpu * 2))

# Load correct modules
module load bowtie2/${bowtie2_version}
module load samtools/${samtools_version}

# Get correct files
get_samples ${wd} ${data_folder}

full_path_fastq="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
fastq=$(basename ${full_path_fastq})

# Check whether script needs to run
if [[ -s "${wd}/data/processed/bowtie2/${sample_id}/${sample_id}_contaminants.sam" ]]; then
  echo "`date` ${sample_id} contaminants already present"
  exit 0
fi

echo "`date` removing contaminants from ${sample_id}"

# Create output dirs
cd "${wd}/data/processed"
mkdir -p "bowtie2/${sample_id}/"

# Run bowtie to remove contaminants
bowtie2 --seedlen=25 \
  --threads ${threads} \
  --time \
  --un-gz "${wd}/data/processed/bowtie2/${sample_id}/${sample_id}_filtered.fastq.gz" \
  -x ${bowtie2_index} \
  -U "${wd}/data/processed/trimgalore/${sample_id}/${fastq}" \
  -S "${wd}/data/processed/bowtie2/${sample_id}/${sample_id}_contaminants.sam"

# Create contaminant QC file
contaminants_type=$(basename $bowtie2_index)
contaminants_file="${wd}/data/processed/bowtie2/${sample_id}/${sample_id}_${contaminants_type}.txt"

# Get total number of reads
tot_reads=$(zcat "${wd}/data/processed/trimgalore/${sample_id}/${fastq}" | echo $((`wc -l`/4)))
echo -e "RiboseQC run for ${sample_id} on `date` \n" >> "${contaminants_file}"
# Print headers to file
printf '\t%s\t%s\t%s\n' "READ_TYPE" "READS" "PERCENTAGE" >> "${contaminants_file}"
# Print total no. of reads
printf '%s\t%s\t%s\t%s\n' $sample_id "Total" $tot_reads "100" >> "${contaminants_file}"

# For each contaminant type, print absolute and relative number of reads
for contaminant_type in tRNA snRNA snoRNA mtDNA rRNA ; do  
  contaminant_reads=`samtools view -@ ${threads} "${wd}/data/processed/bowtie2/${sample_id}/${sample_id}_contaminants.sam" | grep -o "$contaminant_type" | wc -l`
  contaminant_percentage=`awk -vn=248 "BEGIN{print(${contaminant_reads}/${tot_reads}*100)}"`
  printf '%s\t%s\t%s\t%s\n' "${sample_id}" "${contaminant_type}" "${contaminant_reads}" "${contaminant_percentage}" >> "${contaminants_file}"
done

# Count reads that passed filtering
filtered_reads=$(zcat "${wd}/data/processed/bowtie2/${sample_id}/${sample_id}_filtered.fastq.gz" | echo $((`wc -l`/4)))  
filtered_percentage=`awk -vn=248 "BEGIN{print(${filtered_reads}/${tot_reads}*100)}"`
printf '%s\t%s\t%s\t%s\n\n' "${sample_id}" "Passed" "${filtered_reads}" "${filtered_percentage}" >> "${contaminants_file}"

echo "`date` finished ${sample_id}"
