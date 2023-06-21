#!/bin/bash

# Load modules
module load R/${r_version}

# Check whether script needs to run
if [[ -f "${outdir}/customannotation/${merged_gtf_basename}_novel_filtered.gtf" ]]; then
  echo "Merged filtered annotated GTF already present"
  exit 0
fi

# Create output dirs
mkdir -p "${outdir}/customannotation"

echo "`date` Filter and annotate novel GTF"

# Process and filter novel GTF
Rscript "${scriptdir}/filter_annotate.R" \
  "${wd}" \
  "${reference_gtf}" \
  "${outdir}/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}_gffcompare.annotated.gtf" \
  "${refseq_gtf}" \
  "${outdir}/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}_matching.tracking" \
  "${min_occurrence}" \
  "${min_tpm}" \
  "${outdir}/customannotation/${merged_gtf_basename}_novel_filtered.gtf"

echo "`date` Finished GTF filtering"