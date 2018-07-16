library(plyr)
library(doMC)
library(tidyverse)
library(SRAdb)
library(DBI)
library(synapser)

# local
# home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"
# sra_db    <- "/home/aelamb/SRAmetadb.sqlite"

#ec2
home_dir  <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/ubuntu/tmp/"
sra_db    <- "/home/ubuntu/tmp/SRAmetadb.sqlite"

study_id   <- "SRP051688"
upload_id  <- "syn12649849"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

activity_obj <- Activity(
    name = "upload",
    description = "upload raw fastq files from SRA",
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE74246/transfer_fastq_files_to_synapse.R"),
    used = list("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246")
)


if(!file.exists(sra_db)){
    sra_db <- getSRAdbFile()
}

con <- dbConnect(RSQLite::SQLite(), sra_db)


sra_df <- study_id %>%
    listSRAfile(con) %>%
    select(sample, run)

submission_id <- study_id %>% 
    str_c("select * from study where study_accession = '", ., "'") %>% 
    dbGetQuery(con, .) %>% 
    use_series(submission_accession)

query_df <- submission_id %>% 
    str_c(
        "select sample_alias, sample_accession, sample_attribute from sample ",
        "where submission_accession = '", 
        ., 
        "'") %>% 
    dbGetQuery(con, .) %>% 
    left_join(sra_df, by = c("sample_accession" = "sample"))

manifest_df <- query_df %>% 
    set_colnames(c("sample", "SRS_id", "sample_attribute", "SRR_id")) %>% 
    select(sample, SRS_id, SRR_id, sample_attribute) %>% 
    separate(
        sample_attribute,
        sep = " \\|\\| ",
        into = c("cell_type", "day", "cell_description", "assay", "patient"),
        extra = "drop") %>%
    mutate(cell_type = str_remove_all(cell_type, "source_name: ")) %>% 
    mutate(day = str_remove_all(day, "time: ")) %>% 
    mutate(day = str_remove_all(day, " d")) %>% 
    mutate(patient = str_remove_all(patient, "donor: Donor: ")) %>% 
    select(-c(cell_description, assay))
    
transfer_to_synapse <- function(sra){
    getFASTQfile(sra, con)
    fastqs <- list.files(full.names = T) %>% 
        keep(str_detect(., "fastq"))
    walk(fastqs, upload_file_to_synapse, upload_id, activity_obj = activity_obj)
    walk(fastqs, file.remove)
}

file_df <- upload_id %>% 
    get_file_df_from_synapse_dir_id %>% 
    filter(str_detect(file.name, ".fastq.gz$")) %>% 
    mutate(SRR_id = str_match(file.name, "(SRR[0-9]*)_[0-9]{1}.fastq.gz")[,2]) %>% 
    mutate(pair = str_match(file.name, "SRR[0-9]*_([0-9]{1}).fastq.gz")[,2])

manifest_df %>% 
    filter(day == "0") %>% 
    filter(!SRR_id %in% file_df$SRR_id) %>% 
    use_series(SRR_id) %>% 
    walk(transfer_to_synapse)




