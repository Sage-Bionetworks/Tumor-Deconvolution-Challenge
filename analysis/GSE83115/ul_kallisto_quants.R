library(plyr)
library(tidyverse)
library(synapser)
library(doMC)

#ec2
home_dir  <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/ubuntu/tmp/"
kallisto_dir <- "/home/ubuntu/kallisto_cwl/"


manifest_id          <- "syn12177468"
file_view_id         <- "syn12179146"
upload_id            <- "syn12180291"


fasta_file <- "ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
scripts <- list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/dl_fastqs_for_kallisto.R",
                "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/create_kallisto_quants.R",
                "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/ul_kallisto_quants.R")



#ec2
tmp_dir   <- "/home/ubuntu/tmp/"
kallisto_dir <- "/home/ubuntu/kallisto_cwl/"
fastq_file <- "fastq.tsv"
index_file <- "../GRCH38.idx"    

setwd(tmp_dir)
registerDoMC(cores = detectCores() - 1)


file_df <- read_tsv(fastq_file)


ul_quant_files <- function(df){
    abundance_file <- str_c(df$sample_name[[1]], "_kallisto_abundance.tsv")
    activity_obj <- Activity(
        name = "upload",
        description = "upload kallisto quant files",
        executed = list(scripts),
        used = flatten(list(manifest_id, file_view_id, fasta_file, df$id))
    )
    upload_file_to_synapse(abundance_file, upload_id, activity_obj = activity_obj)
}

file_df %>% 
    split(.$sample_name) %>% 
    l_ply(ul_quant_files, .parallel = T)











