###########################################################################
#
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Date: 24-06-2022
#
# Description:
# Script to prepare annotation files for RiboseQC and ORFquant based on
# custom .gtf file created in detect isoform pipeline
#
###########################################################################

# Load libraries ----------------------------------------------------------
message("Loading required libraries ...")
suppressPackageStartupMessages({
library(RiboseQC)
library(rtracklayer)
})

# Set global variables ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

twobit_file <- args[1]
merged_gtf <- args[2]  
annot_dir <- args[3]  
annot_name <- args[4]

# Prepare annotation files ------------------------------------------------
prepare_annotation_files(annotation_directory = annot_dir, 
                         twobit_file = twobit_file, 
                         gtf_file = merged_gtf, 
                         annotation_name = annot_name, 
                         forge_BSgenome = TRUE)

BSgenome_dir <- grep("BSgenome", x = list.dirs(annot_dir, recursive = F), value = T)
BSgenome_package <- basename(BSgenome_dir)


install.packages(BSgenome_dir, 
                character.only = TRUE, 
                repos = NULL, 
                type="source")