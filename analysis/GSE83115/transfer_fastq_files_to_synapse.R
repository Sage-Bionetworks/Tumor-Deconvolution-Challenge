library(plyr)
library(doMC)
library(tidyverse)
library(SRAdb)
library(DBI)
library(synapser)

# local
home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"

#ec2
# home_dir  <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
# tmp_dir   <- "/home/ubuntu/tmp/"

manifest_id <- "syn12177468"
upload_id <- "syn12177447"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

activity_obj <- Activity(
    name = "upload",
    description = "upload raw fastq files from SRA",
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/transfer_fastq_files_to_synapse.R")
)

df <- create_df_from_synapse_id(manifest_id)


transfer_to_synapse <- function(sra){
    getFASTQfile(sra, con)
    fastqs <- list.files(full.names = T) %>% 
        keep(str_detect(., "fastq"))
    walk(fastqs, upload_file_to_synapse, upload_id, activity_obj = activity_obj)
    walk(fastqs, file.remove)
}

walk(df$sample, transfer_to_synapse)



