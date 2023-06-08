######################################################################
#
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
# This script parses merged GTF file as a GRanges object and filters
# and annotates the transcripts based on certain criteria:
#
# - Transcript occurence in sample set
# - Reference gene overlap
# - Novel - ref exon sum difference
# - antisense class i overlap
# - XR / NR overlap annotation
# - Rename transcripts with no ref gene
#
######################################################################

colnames_gtf = c(
  "type",
  "source",
  "gene_id",
  "gene_name",
  "gene_biotype",
  "transcript_id",
  "transcript_name"
)
colnames_stringtie_gtf = c(
  "type",
  "source",
  "gene_id",
  "gene_name",
  "transcript_id",
  "exon_number",
  "ref_gene_id",
  "class_code",
  "cmp_ref"
)

# Load required packages --------------------------------------------------
message(paste(Sys.time(), "Loading required libraries ..."), sep = "\t")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
})

# Global arguments --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

wd = args[1]
gtf_ref_loc = args[2]
gtf_novel_loc = args[3]
gtf_refseq_basename = args[4]
transcript_tracking_loc = args[5]
min_occurence = args[6]
min_tpm = args[7]
outfile = args[8]

# Functions ---------------------------------------------------------------

count_mono_exonics <- function(gtf) {
  # Input:
  # gtf_novel = Granges object containing novel transcripts
  #
  # Finds transcripts with a single exon and removes those from the object
  
  # Remove transcripts that are found by StringTie with only single exon
  t <- as.data.frame(table(gtf$transcript_id))
  colnames(t) <- c("transcript_id",
                   "exon_count")
  
  # Each transcript has a transcript row and exon row
  t$exon_count <- t$exon_count - 1
  gtf_count <- left_join(as.data.frame(gtf),
                         t,
                         by = "transcript_id")
  
  mono_transcripts <- subset(gtf_count,
                             exon_count < 2 &
                               is.na(ref_gene_id) &
                               source == "StringTie")
  
  return(mono_transcripts)
  
}

###########################################################################

calculate_tx_occurence <-
  function(transcript_tracking_loc,
           tracking_file,
           min_occurence = 1,
           min_tpm = 1,
           gtf) {
    # Input:
    # tracking_file = table output from GFFcompare that list occurences of each
    #                 transfrag in each sample
    # gtf_novel = Granges output from GFFcompare GTF
    # min_occurence = integer (default = 1) how many samples should share a transfrag
    # min_tpm = integer (default = 1) how much coverage each sample should have to count as covered
    #
    # Calculate the occurence in the pool of samples, subsets transcripts on
    # minimal occurence parameter. Output occurence file in directory of tracking file.
    # for further filtering.
    
    colnames(tracking_file) <-
      c("transfrag_id", "locus_id", "ref_gene_id", "class_code")
    tracking_file <- subset(tracking_file,
                            class_code == "=")
    
    tracking_genes <- data.frame(samples = tracking_file[, 1:4])
    tracking_genes$rowid <- 1:nrow(tracking_genes)
    
    tracking_samples <-
      data.frame(samples = tracking_file[, 5:length(tracking_file)])
    # tracking_samples$rowid <- 1:nrow(tracking_samples)
    
    tracking_samples_tpm <- sapply(tracking_samples,
                                   function(x)
                                     sapply(strsplit(as.character(x), split = "[|]"), "[", 5))
    
    tracking_samples_tpm <- data.frame(tracking_samples_tpm)
    tracking_samples_tpm <-
      dplyr::mutate_all(tracking_samples_tpm, function(x)
        as.numeric(as.character(x)))
    
    # Count transfrag occurrence
    occurrence_df <-
      data.frame(ifelse(tracking_samples_tpm > min_tpm, 1, 0))
    occurrence_df[is.na(occurrence_df)] <- 0
    occurrence_df <- data.frame(rowSums(occurrence_df))
    colnames(occurrence_df) <- "occurrence"
    sum(occurrence_df > 2)
    
    # Subset transcripts on their occurence
    tracking_genes <- cbind(tracking_genes, occurrence_df)
    tracking_genes <-
      tracking_genes[tracking_genes$occurrence >= min_occurence, ]
    tracking_genes$transcript_id <-
      gsub(".*\\|", "", tracking_genes$samples.ref_gene_id)
    
    # Grab colnames for sample matrix
    sample_ids <- list.files(paste(wd, "data/raw", sep = "/"))
    sample_ids <- unique(gsub("_.*$", "", sample_ids))
    
    # Subset sample matrix and write to tracking file directory
    tracking_file <-
      tracking_file[tracking_file$ref_gene_id %in% tracking_genes$samples.ref_gene_id, ]
    
    colnames(tracking_file) <-
      c("transfrag_id",
        "locus_id",
        "ref_gene_id",
        "class_code",
        sample_ids)
    write.table(
      tracking_file,
      file = gsub(
        "matching.tracking",
        "sample_matrix.txt",
        transcript_tracking_loc
      ),
      quote = F,
      row.names = F,
      sep = "\t"
    )
    
    gtf_keep <- subset(gtf, gtf$transcript_id %in%
                         tracking_genes$transcript_id)
    
    
    return(gtf_keep)
  }

###########################################################################

check_tx_overlap <- function(gtf, gtf_reference) {
  # gtf_novel_processed = input from previous function
  # gtf_reference = ensembl reference GTF Grange
  #
  # This function looks for exonic overlap between query (novel tx exons)
  # and subject (ref tx exons). Novel transcripts that overlap with multiple
  # reference genes are flagged for removal.
  
  # These genes should have no reference overlap
  gtf_no_overlap <- gtf[gtf$class_code %in% c("u", "x", "i", "y"), ]
  no_overlap_tx <- gtf_no_overlap$transcript_id
  
  # Create Granges with reference transcript exons using Ensembl canonical transcript
  gene_tx_refs <-
    gtf_reference[which(gsub(".*-", "", gtf_reference$transcript_name) == "201"), ]$transcript_id
  gene_tx_refs <- unique(gene_tx_refs)
  
  gene_ref_gtf <-
    gtf_reference[which(gtf_reference$type == "exon" &
                          gtf_reference$transcript_id %in% gene_tx_refs), ]
  
  # Grab novel exons
  gtf_novel_exons <-
    GenomicRanges::makeGRangesFromDataFrame(gtf[which(gtf$type == "exon"), ],
                                            keep.extra.columns = T)
  
  gtf_no_overlap_exons <-
    GenomicRanges::makeGRangesFromDataFrame(gtf[which(gtf$type == "exon" &
                                                        gtf$transcript_id %in% no_overlap_tx), ],
                                            keep.extra.columns = T)
  
  # Calculate overlapping exons
  no_overlap <-
    GenomicRanges::findOverlaps(query = gtf_no_overlap_exons,
                                subject = gene_ref_gtf,
                                type = "any")
  if (!(length(no_overlap) == 0)) {
    print("WARNING: Overlap in unannotated genes detected")
    print(length(no_overlap))
  }
  print(length(no_overlap))
  
  overlap <- GenomicRanges::findOverlaps(query = gtf_novel_exons,
                                         subject = gene_ref_gtf,
                                         type = "any")
  overlap_df <- data.frame(gtf_novel_exons[from(overlap),
                                           c("transcript_id",
                                             "type",
                                             "class_code",
                                             "cmp_ref",
                                             "ref_gene_id")],
                           gene_ref_gtf[to(overlap), c("gene_id", "gene_name")])
  
  # Calculate overlapping ref genes
  split_by_tx <- split(overlap_df, overlap_df$transcript_id)
  
  check <- bind_rows(lapply(split_by_tx, function(x) {
    # Count per transcript the number of unique gene IDs with overlapping
    # exons. Label transcripts based on this.
    
    single_check = ifelse(length(unique(x$gene_id)) > 1, "multi_gene", "single_gene")
    df = data.frame(transcript_id = unique(x$transcript_id), single_check)
    
    return(df)
    
  }))
  
  check_no_overlap <- unique(subset(check,
                                    check$single_check == "multi_gene")$transcript_id)
  
  return(check_no_overlap)
}

###########################################################################

filter_i_class <- function(gtf_df, reference_granges) {
  # Input:
  # reference_granges = GRanges object of the reference GTF
  # gtf_novel_df = data frame or data table containing new transcripts
  #
  # Searches for same-strand overlap with known genes for i class
  # transcripts.
  tx_i <-
    gtf_df[which(gtf_df$type == "transcript" &
                   gtf_df$class_code == "i")]$transcript_id
  gtf_i <- gtf_df[which(gtf_df$transcript_id %in% tx_i)]
  i_by_tx <- split(gtf_i,
                   gtf_i$transcript_id)
  
  i_no_pass <- unlist(lapply(i_by_tx, function(x) {
    overlap <-
      GenomicRanges::findOverlaps(
        query = GenomicRanges::makeGRangesFromDataFrame(x,
                                                        keep.extra.columns = T),
        subject = reference_granges,
        type = "any"
      )
    
    if (length(overlap) > 0) {
      return(unique(x$transcript_id))
    }
    
  }))
  
  return(i_no_pass)
  
}

###########################################################################

annotate_overlap <- function(gtf, gtf_refseq_basename, x_name) {
  # Input:
  # gtf_novel = Previous Granges output
  # gtf_refseq_basename = name of custom refseq gtf that holds XR and NR transcripts respectively
  # x_name = character vector of the name of the new column and the name of
  #          the RefSeq GFF that contains the transcripts
  #
  # Annotates transcripts found in gtf_novel with any overlap in the X
  # granges object. Currently used to annotate XR_### and NR_### transcripts
  # as these might not be applicable for neo-antigen detection.
  
  x <-
    rtracklayer::import(paste(gtf_refseq_basename, x_name, "gff", sep = "."))
  
  # Convert refseq seqnames to ensembl seqnames
  x <-
    x[seqnames(x) %in% c(
      "NC_000001.11",
      "NC_000002.12",
      "NC_000003.12",
      "NC_000004.12",
      "NC_000005.10",
      "NC_000006.12",
      "NC_000007.14",
      "NC_000008.11",
      "NC_000009.12",
      "NC_000010.11",
      "NC_000011.10",
      "NC_000012.12",
      "NC_000013.11",
      "NC_000014.9",
      "NC_000015.10",
      "NC_000016.10",
      "NC_000017.11",
      "NC_000018.10",
      "NC_000019.10",
      "NC_000020.11",
      "NC_000021.9",
      "NC_000022.11",
      "NC_000023.11",
      "NC_000024.10"
    )]
  seqlevels(x) <- as.character(unique(seqnames(x)))
  x <- GenomeInfoDb::renameSeqlevels(x, c(1:22, "X", "Y"))
  x <- subset(x, x$type == "exon")
  
  gtf_novel_gr <-
    GenomicRanges::makeGRangesFromDataFrame(gtf, keep.extra.columns = T)
  
  gtf_novel_exons <-
    subset(gtf_novel_gr, gtf_novel_gr$type == "exon")
  
  x_overlap <-
    GenomicRanges::findOverlaps(query = x,
                                subject = gtf_novel_exons,
                                type = "any")
  
  novel_x_hits <-
    as.data.frame(unique(gtf_novel_exons[subjectHits(x_overlap)]))
  novel_x_hits <- unique(novel_x_hits[, "transcript_id"])
  elementMetadata(gtf_novel_gr)[[paste0(x_name, "_overlap")]] <-
    ifelse(gtf_novel_gr$transcript_id %in% novel_x_hits,
           paste0(x_name, "_hit"),
           "none")
  return(gtf_novel_gr)
}

###########################################################################

rename_stringtie_transcripts <- function(gtf_novel_df) {
  # Input:
  # gtf_novel_df = data frame or data table containing new transcripts
  #
  # Renames genes without a reference gene in a similar fashion as
  # we name our ORFs using chromosome, start and stop position. Linked
  # transcripts will have the same gene name
  
  gtf_by_gene <- split(gtf_novel_df, gtf_novel_df$gene_id)
  
  gtf_renamed <-
    suppressWarnings(bind_rows(lapply(gtf_by_gene, function(x) {
      if (is.na(x$ref_gene_id)[1]) {
        chr <- unique(x$seqnames)
        start <- min(x$start)
        end <- max(x$end)
        x$gene_name <- paste0(chr, ":", start, "-", end)
      }
      
      return(x)
      
    })))
}

# Load required files -----------------------------------------------------
message(paste(Sys.time(), "Importing GTF files ..."), sep = "\t")
gtf_reference <-
  rtracklayer::import.gff(gtf_ref_loc, colnames = colnames_gtf)
gtf_ref_df <- as.data.frame(gtf_reference)
gtf_novel_file <-
  rtracklayer::import.gff(gtf_novel_loc, colnames = colnames_stringtie_gtf)
tracking_file <-
  read.table(file = transcript_tracking_loc, header = F)

# Filtering and annotation ------------------------------------------------
message(paste(Sys.time(), "Removing transcripts on scaffolds ..."), sep = "\t")
gtf_novel <-
  gtf_novel_file[seqnames(gtf_novel_file) %in% c(1:22, "X", "Y"), ]

message(paste(Sys.time(), "Removing unstranded transcripts ..."), sep = "\t")
gtf_novel <- gtf_novel[strand(gtf_novel) != "*", ]

message(paste(Sys.time(), "Removing mono-exonic transcripts ..."), sep = "\t")
mono_exonic_novel <- count_mono_exonics(gtf = gtf_novel)
gtf_novel <-
  gtf_novel[!(gtf_novel$transcript_id %in% mono_exonic_novel$transcript_id), ]

# Flagging RefSeq transcripts ---------------------------------------------
message(paste(Sys.time(), "Flagging XR transcript overlap ..."), sep = "\t")
gtf_novel <- annotate_overlap(gtf = gtf_novel,
                              gtf_refseq_basename = gtf_refseq_basename,
                              x_name = "xr")

message(paste(Sys.time(), "Flagging NR transcript overlap ..."), sep = "\t")
gtf_novel <- annotate_overlap(gtf = gtf_novel,
                              gtf_refseq_basename = gtf_refseq_basename,
                              x_name = "nr")

# Add biotype to custom annotation based on reference ID ------------------
message(paste(Sys.time(), "Adding biotype to StringTie transcripts ..."),
        sep = "\t")
gtf_novel$gene_biotype <-
  gtf_reference$gene_biotype[match(gtf_novel$ref_gene_id,
                                   gtf_reference$gene_id)]

# Remove unwanted transcripts ---------------------------------------------
transcripts_discard <-
  gtf_novel[which(
    gtf_novel$type == "transcript" &
      gtf_novel$class_code %in% c("=", "c", "j", "m", "n", "e", "r", "s")
  ), ]
gtf_novel_keep <-
  gtf_novel[which(!(
    gtf_novel$transcript_id %in% transcripts_discard$transcript_id
  )), ]

# Transcript occurence ----------------------------------------------------
message(paste(Sys.time(), "Counting transcript occurence in samples ..."),
        sep = "\t")
gtf_novel_keep <-
  calculate_tx_occurence(
    transcript_tracking_loc = transcript_tracking_loc,
    tracking_file = tracking_file,
    min_occurence = min_occurence,
    gtf = gtf_novel_keep
  )

gtf_novel_df <- as.data.frame(gtf_novel_keep)

# Annotate gene names, reference IDs and biotypes of exons ----------------
message(paste(Sys.time(), "Annotate gene names, reference IDs and biotypes ..."),
        sep = "\t")
gtf_by_tx <- split(gtf_novel_df, gtf_novel_df$transcript_id)

gtf_novel_processed <- bind_rows(lapply(gtf_by_tx, function(x) {
  # Find existing ensembl gene IDs
  # gene_name and ref_gene_id are annotated in the first row of each transcript
  
  if (!is.na(x$gene_name)[1]) {
    x$gene_name <- x$gene_name[1]
  }
  if (!is.na(x$ref_gene_id[1])) {
    x$ref_gene_id <- x$ref_gene_id[1]
  }
  if (!is.na(x$class_code[1])) {
    x$class_code <- x$class_code[1]
  }
  x$gene_biotype <- x$gene_biotype[1]
  
  known_ref_ids <-
    length(unique(x$ref_gene_id[grep("ENS", x$ref_gene_id)]))
  
  if (known_ref_ids == 1) {
    # fill the missing gene_id, gene_name and gene_biotype into novel transcripts
    x$gene_id      <-
      unique(x$ref_gene_id[grep("ENS", x$ref_gene_id)])
    x$gene_biotype <-
      unique(x$gene_biotype[grep("ENS", x$ref_gene_id)])
    
  }
  
  else if (known_ref_ids > 1) {
    # if reference is know then replace the gene ids
    x$gene_id[grep("ENS", x$ref_gene_id)] = x$ref_gene_id[grep("ENS", x$ref_gene_id)]
    
    # which transcripts have to be removed
    coord_remove <- grep("ENS", x$ref_gene_id, invert = T)
    
    # remove the unwanted transcripts
    if (any(coord_remove)) {
      x <- x[-coord_remove, ]
    }
    
  }
  
  else if (known_ref_ids == 0) {
    x$gene_name    <- x$gene_id
    x$gene_biotype <- "stringtie"
  }
  
  return(x)
  
}))

# Check ref tx overlap ----------------------------------------------------
message(paste(Sys.time(), "Checking reference transcript overlap ..."),
        sep = "\t")
ref_overlap_txs <- check_tx_overlap(gtf = gtf_novel_processed,
                                    gtf_reference = gtf_reference)

# Annotate same-stranded i class ------------------------------------------
message(paste(Sys.time(), "Remove same sense i class transcripts ..."),
        sep = "\t")
same_strand_i_txs <- filter_i_class(gtf_df = gtf_novel_keep,
                                    reference_granges = gtf_reference)

# Filter for final transcript subset --------------------------------------
gtf_novel_filtered <- subset(
  gtf_novel_processed,!(gtf_novel_processed$transcript_id %in% same_strand_i_txs) &
    !(gtf_novel_processed$transcript_id %in% ref_overlap_txs)
)

gr_novel_filtered <-
  GenomicRanges::makeGRangesFromDataFrame(gtf_novel_filtered, keep.extra.columns = T)

export(object = gr_novel_filtered,
       con = outfile)

# Rename new stringtie genes ----------------------------------------------
message(paste(Sys.time(), "Rename novel genes ..."), sep = "\t")
gtf_novel_processed_final <-
  rename_stringtie_transcripts(gtf_novel_df = gtf_novel_filtered)
gtf_novel_GR <-
  GenomicRanges::makeGRangesFromDataFrame(gtf_novel_processed_final, keep.extra.columns = T)

# Merge novel GTF with hg38 GTF -------------------------------------------
gtf_novel_merged <- c(gtf_reference, gtf_novel_GR)

# Export GTF --------------------------------------------------------------
message(paste(Sys.time(), "Exporting custom gtf ... "), sep = "\t")
rtracklayer::export.gff(object = gtf_novel_merged, con = outfile)
message(paste(Sys.time(), "Exporting custom gtf ... Done!"), sep = "\t")

# Plot class code distribution --------------------------------------------
message(paste(Sys.time(), "Plotting transcript classes ... "), sep = "\t")
class_data <- gtf_unique_df[gtf_unique_df$type == "transcript" &
                              !(is.na(gtf_unique_df$class_code)), c("class_code", "transcript_id")]

ggplot(class_data) +
  geom_bar(aes(x = class_code, y = stat(count), fill = class_code)) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  device = "pdf",
  filename = "tx_class_distribution_all.pdf",
  width = 8,
  height = 8,
  path = paste(wd, "data/processed", sep = "/")
)

ggplot(class_data[class_data$class_code %in% c("u", "x", "y"), ]) +
  geom_bar(aes(x = class_code, y = stat(count), fill = class_code)) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  device = "pdf",
  filename = "tx_class_distribution_novel.pdf",
  width = 8,
  height = 8,
  path = paste(wd, "data/processed", sep = "/")
)

message(paste(Sys.time(), "Finished filter & annotate steps "), sep = "\t")