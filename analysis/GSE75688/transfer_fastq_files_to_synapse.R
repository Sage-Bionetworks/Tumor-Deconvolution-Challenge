library(SRAdb)
library(DBI)
library(synapser)

home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE75688/"
sra_id    <- "SRP102688"
upload_id <- "syn11969192"
sra_db    <- "/home/aelamb/SRAmetadb.sqlite"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

activity_obj <- Activity(
    name = "upload",
    description = "upload raw fastq files from SRA",
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE75688/transfer_fastq_files_to_synapse.R")
)

if(!file.exists(sra_db)){
    sra_db <- getSRAdbFile()
}

con <- dbConnect(RSQLite::SQLite(), sra_db)
df <- listSRAfile(sra_id, con)

transfer_to_synapse <- function(){
    getFASTQfile(df$sample[[1]])
    fastq <- list.files(full.names = T) %>% 
        keep(str_detect(., "fastq"))
    upload_file_to_synapse(fastq, upload_id, activity_obj = activity_obj)
    file.remove(fastq)
}

walk(df$sample, transfer_to_synapse)
