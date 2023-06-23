###########################################################################
#
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Date: 04-06-2021
#
# Description:
# R script to count ORF occurrence across samples, taking the ORFs from 
# the merged analyses as reference input, and looping through the individual 
# samples to find matching ORFs. Exports a .csv file containing some stats 
# (the specific stats to display can be customized) and the number of samples
# with matching ORFs.
#
###########################################################################

# Load libraries ----------------------------------------------------------
message("Loading required libraries ...")
suppressPackageStartupMessages({
library(ORFquant)
library(GenomicRanges)
})

# Define general variables ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
wd = args[1]
ORFquant_results = args[2]
pool_id = args[3]

# Define functions --------------------------------------------------------

load_pooled_results <- function(ORFquant_results) {

    message("Loading pooled ORFquant results ...")
    load(paste0(ORFquant_results))

    pooled_gen_regions <- ORFquant_results$ORFs_gen
    pooled_stats <- ORFquant_results$ORFs_tx

    pooled_list <- split(pooled_gen_regions, names(pooled_gen_regions))
    mcols(pooled_list) <- pooled_stats[names(pooled_list)]

    return(pooled_list)
}

load_sample_results <- function(workdir) {
    
    # Find all ORFquant_results files 
    # (excluding pooled results)

    message("Loading sample ORFquant results ...")

    fnames <- grep(list.files(pattern = "_final_ORFquant_results$", recursive = TRUE, path = workdir, full.names = T), 
                   pattern = pool_id, invert = T, value = T)
    sample_ids <- grep(list.dirs(paste0(workdir, "/analysis"), recursive = F, full.names = F), 
                       pattern = pool_id, invert = T, value = T)
    sample_df <- cbind(fnames, sample_ids)

    return(sample_df)
}

create_sample_granges <- function(fnames, sample_ids) {

    # Load all ORFquant_results files into a GrangesList and a list of dataframes

    message("Loading per-sample ORFquant_results files ...")

    all_ORFquant_results <- GRangesList()

    for (i in seq_along(fnames)) {

        message(paste("Loading", fnames[i], "..."))
        fname = fnames[i]
        sample_id = sample_ids[i]
  
        load(fname)
        gen_regions <- ORFquant_results$ORFs_gen
  
        all_ORFquant_results[[sample_id]] <- gen_regions
        }

    message("Loading per-sample ORFquant_results files ... Done!")

    return(all_ORFquant_results)
}

match_ORFs <- function(sample_ids,all_ORFquant_results, pooled_list) {

    # Create 'master data.frame'
    pooled_df <- data.frame("ORF_id" = names(pooled_list),
                            chr = unique(data.frame(seqnames(pooled_list)))$value,
                            start = min(start(pooled_list)), 
                            end = max(end(pooled_list)), 
                            exon_width = sum(width(pooled_list)),
                            pval = mcols(pooled_list)$pval, 
                            count_exact_match = 0, 
                            count_same_stop = 0,
                            ORFs_pM = mcols(pooled_list)$ORFs_pM,
                            Protein = mcols(pooled_list)$Protein, 
                            biotype_tx = mcols(pooled_list)$transcript_biotype, 
                            biotype_gene = mcols(pooled_list)$gene_biotype,
                            category_tx = mcols(pooled_list)$ORF_category_Tx_compatible,
                            category_gen = mcols(pooled_list)$ORF_category_Gen,
                            gene_id = mcols(pooled_list)$gene_id,
                            gene_name = mcols(pooled_list)$gene_name)

    for (i in seq_along(sample_ids)) {
  
        message(paste0("Matching reference ORFs with ORFs in sample ", sample_ids[i]))

        # Grab ORFquant_results for subject_name and split into list by transcript ID
        subject_sample <- all_ORFquant_results[[sample_ids[i]]]
        subject_list <- split(subject_sample, names(subject_sample))
        
        # Find exact hits
        exact_hits <- findOverlaps(query = pooled_list, subject = subject_list, type = "equal")
        pooled_exact_hits <- pooled_list[queryHits(exact_hits)]
        subject_exact_hits <- subject_list[subjectHits(exact_hits)]
        
        # Calculate overlap
        exact_overlaps <- GenomicRanges::intersect(pooled_exact_hits, subject_exact_hits)
        percent_overlap_exact <- sapply(sum(width(exact_overlaps)) / sum(width(pooled_exact_hits)), max)
        
        # Grab matching ORFs with > 90% overlap
        matching_ORFs_exact <- names(percent_overlap_exact[which(percent_overlap_exact > 0.9)])
        
        pooled_df[which(pooled_df$ORF_id %in% matching_ORFs_exact), ]$count_exact_match <- 
            pooled_df[which(pooled_df$ORF_id %in% matching_ORFs_exact), ]$count_exact_match + 1
        

        # Find hits with same stop
        end_hits <- findOverlaps(query = pooled_list, subject = subject_list, type = "end")
        pooled_end_hits <- pooled_list[queryHits(end_hits)]
        subject_end_hits <- subject_list[subjectHits(end_hits)]
        
        # Calculate overlap
        end_overlaps <- GenomicRanges::intersect(pooled_end_hits, subject_end_hits)
        percent_overlap_end <- sapply(sum(width(end_overlaps)) / sum(width(pooled_end_hits)), max)
        
        # Grab matching ORFs with > 90% overlap
        matching_ORFs_end <- names(percent_overlap_end[which(percent_overlap_end > 0.9)])
        
        pooled_df[which(pooled_df$ORF_id %in% matching_ORFs_end), ]$count_same_stop <-
            pooled_df[which(pooled_df$ORF_id %in% matching_ORFs_end), ]$count_same_stop + 1
    }

    return(pooled_df)

}

# Load data ---------------------------------------------------------------

pooled_list <- load_pooled_results(ORFquant_results = ORFquant_results)

sample_df <- load_sample_results(workdir = wd)

sample_granges <- create_sample_granges(fnames = sample_df[,1],
                                        sample_ids = sample_df[,2])

message("Loading data ... Done!")

# Count ORF occurences ----------------------------------------------------

pooled_df <- match_ORFs(sample_ids = sample_df[,2],
                        all_ORFquant_results = sample_granges,
                        pooled_list = pooled_list)

message("Matching reference ORFs with ORFs in samples ... Done!")

# Export results to CSV ---------------------------------------------------

message("Exporting results to csv file ... ")

write.csv2(x = pooled_df, file = paste0(wd, "/analysis/",pool_id,"/",pool_id,"_ORFcounts.csv"), quote = F)

message("Exporting results to csv file ... Done!")