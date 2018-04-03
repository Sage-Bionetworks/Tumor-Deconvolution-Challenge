library(SRAdb)
library(DBI)
library(synapser)

# local
home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE72056/"
sra_db    <- "/home/aelamb/SRAmetadb.sqlite"

# ec2
# home_dir  <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
# tmp_dir   <- "/home/ubuntu/tmp/"
# sra_db    <- "/home/ubuntu/tmp/SRAmetadb.sqlite"


sra_id    <- "?"
upload_id <- "syn12063033"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

activity_obj <- Activity(
    name = "upload",
    description = "upload raw fastq files from SRA",
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE72056/transfer_fastq_files_to_synapse.R")
)

if(!file.exists(sra_db)){
    sra_db <- getSRAdbFile()
}

con <- dbConnect(RSQLite::SQLite(), sra_db)
df <- listSRAfile(sra_id, con)

transfer_to_synapse <- function(sra){
    getFASTQfile(sra, con)
    fastqs <- list.files(full.names = T) %>% 
        keep(str_detect(., "fastq"))
    walk(fastqs, upload_file_to_synapse, upload_id, activity_obj = activity_obj)
    walk(fastqs, file.remove)
}

walk(df$sample, transfer_to_synapse)
