#!/bin/bash

set -uo pipefail

# Get sample IDs
mapfile -t sample_ids < sample_ids.txt

# Get gtfmergefile
gtfmergefile="${outdir}/stringtie/gtfmergefile.txt"

# Load modules
module load stringtie/${stringtie_version}
module load gffcompare/${gffcompare_version}

# Check whether script needs to run
if [[ -f "${outdir}/stringtie/${merged_gtf_basename}/${merged_gtf_basename}.gtf" ]]; then
  echo "Merged GTF already present"
  exit 0
fi

# Create output dirs
mkdir -p "${outdir}/gffcompare/${merged_gtf_basename}/"
mkdir -p "${outdir}/stringtie/${merged_gtf_basename}/"
cd "${outdir}/stringtie/${merged_gtf_basename}/"

echo "`date` running Stringtie --merge"

# Run stringtie merge
stringtie --version
stringtie ${gtfmergefile} \
  --merge \
  -G ${reference_gtf} \
  -f 0.05 \
  -m 50 \
  -T 5 \
  -o "${outdir}/stringtie/${merged_gtf_basename}/${merged_gtf_basename}.gtf" \
  -p 4

echo "`date` finished Stringtie Merge"

echo "`date` running GFFcompare for annotation"

# Create output dirs
cd "${outdir}/gffcompare/${merged_gtf_basename}/"

# Run GFFcompare to annotate novel assembly .GTF
gffcompare \
  -V \
  -r ${reference_gtf} \
  -s ${masked_fasta} \
  -o "${merged_gtf_basename}_gffcompare" \
  "${outdir}/stringtie/${merged_gtf_basename}/${merged_gtf_basename}.gtf"

echo "`date` running GFFcompare for transcript occurence"

# Run GFFcompare to generate .tracking file

for i in ${!sample_ids[@]}; do
  gtf_list[i]="${outdir}/stringtie/${sample_ids[i]}/${sample_ids[i]}.gtf"
done

gffcompare \
  -V \
  -r "${outdir}/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}_gffcompare.annotated.gtf" \
  -o "${merged_gtf_basename}_matching" \
  ${gtf_list[@]}

echo "`date` finished running GFFcompare"