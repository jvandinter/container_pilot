#!/bin/bash

# Run parameters
export wd=`pwd`
export project_data_folder="${wd}/data/raw/"
export outdir="${wd}/data/processed"
export simul_array_runs=24
export species="Homo_sapiens"
export genome_version="GRCh38"
export annot="102"
export paired_end="true"
export projectname="NB_tumoroids"
export merged_gtf_basename="NB_tumoroids"
export create_annotation="TRUE"

# Transcript filtering parameters
export min_occurrence=3
export min_tpm=1 

# Set paths
export resource_dir="/hpc/pmc_vanheesch/shared_resources/"
export scriptdir="${wd}/scripts/"
export data_folder="/hpc/pmc_vanheesch/data"  # Data folder containing all of our sequencing data

# Set reference files
export star_index_basedir="${resource_dir}/GENOMES/Homo_sapiens.GRCh38/102/STAR/2.7.8a"
export reference_gtf="${resource_dir}/GENOMES/Homo_sapiens.GRCh38/102/annotation/Homo_sapiens.GRCh38.102.gtf"
#export reference_gtf="/hpc/pmc_vanheesch/projects/Damon/Neoantigens_pediatric_cancers/neuroblastoma/20220913_NB_Tumoroids_RNAseq_Entinostat_691B_Annelisa_rerun/data/processed/customannotation/691B_merged_novel_filtered.gtf"
export refseq_gtf="${resource_dir}/GENOMES/Homo_sapiens.GRCh38/102/annotation/Homo_sapiens.GRCh38.p13"
export reference_genome="/${resource_dir}/GENOMES/Homo_sapiens.GRCh38/102/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
export masked_fasta="${resource_dir}/GENOMES/Homo_sapiens.GRCh38/102/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
export twobit="${resource_dir}/GENOMES/Homo_sapiens.GRCh38/102/Homo_sapiens.GRCh38.dna.primary_assembly.2bit"
export kallisto_index="${resource_dir}/GENOMES/Homo_sapiens.GRCh38/102/kallisto/0.44/kallisto_index"

# Module versions
export cutadapt_version=3.4
export fastqc_version=0.11.9
export trimgalore_version=0.6.6
export star_version=2.7.8a
export samtools_version=1.12
export python_version=3.6.1
export subread_version=2.0.2
export stringtie_version=2.1.5
export gffcompare_version=0.12.2
export r_version=4.1.2
export gffread_version=0.12.6
export salmon_version=1.8.0
export multiqc_version=1.11


###########################
### Resource allocation ###
###########################
low_mem=8G
medium_mem=48G
high_mem=200G
low_cpu=1
medium_cpu=6
high_cpu=16