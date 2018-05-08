# created manifest for getting raw fastq files form geo

library(plyr)
library(doMC)
library(tidyverse)
library(SRAdb)
library(DBI)
library(synapser)

# local
home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"
sra_db    <- "/home/aelamb/SRAmetadb.sqlite"

#ec2
# home_dir  <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
# tmp_dir   <- "/home/ubuntu/tmp/"
# sra_db    <- "/home/ubuntu/tmp/SRAmetadb.sqlite"


sra_id    <- "SRP076277"
upload_id <- "syn12177447"

immune_cell_df <- data_frame(
    "id" = c(
        "SRX1830409",
        "SRX1830408",
        "SRX1830407",
        "SRX1830406",
        "SRX1830395",
        "SRX1830394",
        "SRX1830393",
        "SRX1830392"),
    "sample_name" = c(
        "CD8+_T_cells_COH",
        "CD4+_T_cells_COH",
        "Monocytes_COH",
        "B_cells_COH",
        "CD8+_T_cells_BGI",
        "CD4+_T_cells_BGI",
        "Monocytes_BGI",
        "B_cells_BGI")
)


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

activity_obj <- Activity(
    name = "upload",
    description = "upload raw fastq files from SRA",
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/create_fastq_manifest.R")
)

if(!file.exists(sra_db)){
    sra_db <- getSRAdbFile()
}

con <- dbConnect(RSQLite::SQLite(), sra_db)
df <- listSRAfile(sra_id, con) %>% 
    inner_join(immune_cell_df, by = c("experiment" = "id"))

write_tsv(df, "fastq_id.tsv")
upload_file_to_synapse("fastq_id.tsv", upload_id, activity_obj = activity_obj)