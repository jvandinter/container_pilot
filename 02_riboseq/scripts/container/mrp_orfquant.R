###########################################################################
#
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Date: 04-06-2021
#
# Description:
# Script to run ORFquant for ORF detection and quantification
# on sample-to-sample basis
#
###########################################################################

# Load libraries ----------------------------------------------------------
message("Loading required libraries ...")
suppressPackageStartupMessages({
  library(ORFquant)
  library(RiboseQC)
  library(rmarkdown)
})

# Get variables from input ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
wd = args[1]
for_ORFquant_file = args[2]
sample_id = args[3]
annot_name = args[4]
rannot = args[5]
cpu = args[6]
pandoc_dir = args[7]
resource_dir = args[8]
annotation_package = args[9]

# Define variables --------------------------------------------------------
html_report = paste0(wd,"/data/processed/ORFquant/",sample_id,"/",sample_id,"_report.html")
ORFquant_out = paste0(wd,"/data/processed/ORFquant/",sample_id,"/",sample_id)
ORFquant_output_file = paste0(wd,"/data/processed/ORFquant/",sample_id,"/",sample_id,"_final_ORFquant_results")
ORFquant_plot_data = paste0(ORFquant_output_file, "_plots/",sample_id,"_ORFquant_plots_RData")
find_pandoc(dir = pandoc_dir)

# Define functions --------------------------------------------------------
ORFquant_analysis <- function(for_ORFquant_file, ORFquant_output_file, ORFquant_plot_data, rannot, cpu, html_report, sample_id) {

    print(for_ORFquant_file)

    if(!file.exists(ORFquant_output_file)) {

        run_ORFquant(for_ORFquant_file = for_ORFquant_file,
                     annotation_file = rannot,
                     n_cores = cpu,
                     prefix = ORFquant_out,
                     gene_name = NA,
                     gene_id = NA,
                     genomic_region = NA,
                     write_temp_files = T,
                     write_GTF_file = T,
                     write_protein_fasta = T,
                     interactive = T,
                     stn.orf_find.all_starts = T,
                     stn.orf_find.nostarts = F,
                     stn.orf_find.start_sel_cutoff = NA,
                     stn.orf_find.start_sel_cutoff_ave = 0.5,
                     stn.orf_find.cutoff_fr_ave = 0.5,
                     stn.orf_quant.cutoff_cums = NA,
                     stn.orf_quant.cutoff_pct = 2,
                     stn.orf_quant.cutoff_P_sites = NA,
                     unique_reads_only = F,
                     canonical_start_only = T,
                     stn.orf_quant.scaling = "total_Psites"
                     )
    } else {
        print("ORFquant file already exists")
    }

    plot_ORFquant_results(for_ORFquant_file = for_ORFquant_file,
        ORFquant_output_file = ORFquant_output_file,
        annotation_file = rannot
        )

    create_ORFquant_html_report(input_files = ORFquant_plot_data,
        input_sample_name = sample_id,
        output_file= html_report
        )
}

# Run script --------------------------------------------------------------

if(!require(basename(annotation_package), character.only = T)) {
  install.packages(annotation_package,
                   repos=NULL, type="source")
}

ORFquant_analysis(for_ORFquant_file = for_ORFquant_file,
                  ORFquant_output_file = ORFquant_output_file,
                  ORFquant_plot_data = ORFquant_plot_data,
                  rannot = rannot,
                  cpu = cpu,
                  html_report = html_report,
                  sample_id = sample_id)
