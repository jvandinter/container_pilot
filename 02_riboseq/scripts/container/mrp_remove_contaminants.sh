#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 03-06-2021
#
######################################################################

# Load parameters from main script
threads=$((SLURM_CPUS_PER_TASK * 2))

# Load files
mapfile -t r1_files < ${project_folder}/documentation/r1_files.txt
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Set names
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
r1_file="${r1_files[$((SLURM_ARRAY_TASK_ID-1))]}"
r1_filename=$(basename ${r1_file})

# Check whether script needs to run
if [[ -s "${outdir}/bowtie2/${sample_id}/${sample_id}_contaminants.sam" ]]; then
  echo "`date` ${sample_id} contaminants already present"
  exit 0
fi

echo "`date` removing contaminants from ${sample_id}"

# Create output dirs
cd "${outdir}"
mkdir -p "bowtie2/${sample_id}/"

# Run bowtie to remove contaminants
apptainer exec -B /hpc:/hpc ${container_dir}/bowtie2-2.4.2.sif bowtie2 \
  --seedlen=25 \
  --threads ${threads} \
  --time \
  --un-gz "${outdir}/bowtie2/${sample_id}/${sample_id}_filtered.fastq.gz" \
  -x ${bowtie2_index} \
  -U "${outdir}/trimgalore/${sample_id}/${r1_filename}" \
  -S "${outdir}/bowtie2/${sample_id}/${sample_id}_contaminants.sam"

# Create contaminant QC file
contaminants_type=$(basename $bowtie2_index)
contaminants_file="${outdir}/bowtie2/${sample_id}/${sample_id}_${contaminants_type}.txt"

# Get total number of reads
tot_reads=$(zcat "${outdir}/trimgalore/${sample_id}/${r1_filename}" | echo $((`wc -l`/4)))
echo -e "RiboseQC run for ${sample_id} on `date` \n" >> "${contaminants_file}"
# Print headers to file
printf '\t%s\t%s\t%s\n' "READ_TYPE" "READS" "PERCENTAGE" >> "${contaminants_file}"
# Print total no. of reads
printf '%s\t%s\t%s\t%s\n' $sample_id "Total" $tot_reads "100" >> "${contaminants_file}"

# For each contaminant type, print absolute and relative number of reads
for contaminant_type in tRNA snRNA snoRNA mtDNA rRNA ; do  
  contaminant_reads=`apptainer exec -B /hpc:/hpc ${container_dir}/samtools-1.12.sif samtools view -@ ${threads} "${outdir}/bowtie2/${sample_id}/${sample_id}_contaminants.sam" | grep -o "$contaminant_type" | wc -l`
  contaminant_percentage=`awk -vn=248 "BEGIN{print(${contaminant_reads}/${tot_reads}*100)}"`
  printf '%s\t%s\t%s\t%s\n' "${sample_id}" "${contaminant_type}" "${contaminant_reads}" "${contaminant_percentage}" >> "${contaminants_file}"
done

# Count reads that passed filtering
filtered_reads=$(zcat "${outdir}/bowtie2/${sample_id}/${sample_id}_filtered.fastq.gz" | echo $((`wc -l`/4)))  
filtered_percentage=`awk -vn=248 "BEGIN{print(${filtered_reads}/${tot_reads}*100)}"`
printf '%s\t%s\t%s\t%s\n\n' "${sample_id}" "Passed" "${filtered_reads}" "${filtered_percentage}" >> "${contaminants_file}"

echo "`date` finished ${sample_id}"
