library(SRAdb)
library(DBI)
library(synapser)

home_dir  <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/tmp/"
sra_id    <- "SRP048971"
upload_id <- "syn11953213"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

activity_obj <- Activity(
    name = "upload",
    description = "upload raw fastq files from SRA",
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE62408/transfer_fastq_files_to_synapse.R")
)

srafile = getSRAdbFile()
#srafile <- "/home/aelamb/SRAmetadb.sqlite"
con <- dbConnect(RSQLite::SQLite(), srafile)
df <- listSRAfile(sra_id, con)
getFASTQfile(df$sample, con)
files <- list.files(full.names = T)
walk(files, upload_file_to_synapse, upload_id, activity_obj = activity_obj)
