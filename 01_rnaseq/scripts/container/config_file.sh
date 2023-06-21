#!/bin/bash

# Run parameters
export wd=""
export project_data_folder="${wd}/data/raw/"
export outdir="${wd}/data/processed"
export simul_array_runs=24
export paired_end="true"
export projectname=""
export merged_gtf_basename=""
export create_annotation="TRUE"

# Reference parameters
export species="Homo_sapiens"
export genome_version="GRCh38"
export annot_version="102"

# Transcript filtering parameters
export min_occurrence=3
export min_tpm=1 

# Set paths
export resource_dir="/hpc/pmc_vanheesch/shared_resources"
export scriptdir="${wd}/container_pilot/RNAseq_container" # Folder containing all scripts for the pipeline
export container_dir="/hpc/local/Rocky8/pmc_vanheesch/singularity_images" # Folder containing all of our container versions
export data_folder="/hpc/pmc_vanheesch/data"  # Data folder containing all of our sequencing data

# Set reference files
export star_index_basedir="${resource_dir}/GENOMES/${species}.${genome_version}/${annot_version}/STAR/2.7.8a"
export reference_gtf="${resource_dir}/GENOMES/${species}.${genome_version}/${annot_version}/annotation/${species}.${genome_version}.${annot_version}.gtf"
export refseq_gtf="${resource_dir}/GENOMES/${species}.${genome_version}/${annot_version}/annotation/${species}.${genome_version}.p13"
export reference_genome="/${resource_dir}/GENOMES/${species}.${genome_version}/${annot_version}/${species}.${genome_version}.dna.primary_assembly.fa"
export masked_fasta="${resource_dir}/GENOMES/${species}.${genome_version}/${annot_version}/${species}.${genome_version}.dna_sm.primary_assembly.fa"
export twobit="${resource_dir}/GENOMES/${species}.${genome_version}/${annot_version}/${species}.${genome_version}.dna.primary_assembly.2bit"
export kallisto_index="${resource_dir}/GENOMES/${species}.${genome_version}/${annot_version}/kallisto/0.44/kallisto_index"
