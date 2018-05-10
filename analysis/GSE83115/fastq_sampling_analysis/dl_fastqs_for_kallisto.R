library(plyr)
library(tidyverse)
library(synapser)
library(doMC)

# local
# home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# work_dir   <- "/home/aelamb/"

# ec2
home_dir     <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
work_dir     <- "/home/ubuntu/"

manifest_id          <- "syn12177468"
file_view_id         <- "syn12179146"
download_id_raw      <- "syn12177447"
download_id_sampled  <- "syn12177448"
index_id             <- "syn12212935"

cell_lines <- c("Monocytes_BGI", "CD4+_T_cells_BGI", "CD8+_T_cells_BGI")


setwd(home_dir)
source("scripts/utils.R")
setwd(work_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

manifest_df <- manifest_id %>% 
    create_df_from_synapse_id %>% 
    filter(sample_name %in% cell_lines)

file_view_df <- synTableQuery("select id, name, parentId from syn12179146") %>% 
    as.data.frame %>% 
    filter(str_detect(name, "fastq.gz$")) %>% 
    select(id, name, parentId)

raw_file_df <- file_view_df %>% 
    filter(parentId == download_id_raw) %>% 
    mutate(run = str_sub(name, end = 10)) %>% 
    inner_join(manifest_df) %>% 
    select(name, id, sample_name) %>% 
    mutate(pair = str_match(name, "^[:print:]+_([1,2]).fastq.gz$")[,2])


sampled_file_df <- file_view_df %>% 
    filter(parentId == download_id_sampled) %>% 
    mutate(pair = str_match(name, "^[:print:]+_p([1,2]).fastq.gz$")[,2]) %>% 
    mutate(sample_name = str_match(name, "^([:print:]+)_p[1,2].fastq.gz$")[,2])

file_df <- list(raw_file_df, sampled_file_df) %>% 
    bind_rows %>% 
    select(-parentId)

file_df <- list(raw_file_df, sampled_file_df) %>% 
    bind_rows %>% 
    select(-parentId) 

paths <- llply(file_df$id, download_from_synapse, .parallel = T)

file_df$paths <- unlist(paths)

write_tsv(file_df, "fastq.tsv")

    













