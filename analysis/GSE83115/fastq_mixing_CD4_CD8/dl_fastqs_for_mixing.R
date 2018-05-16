library(plyr)
library(doMC)
library(tidyverse)
library(magrittr)
library(synapser)

# local
# home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"

#ec2
home_dir  <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/ubuntu/tmp/"

manifest_id  <- "syn12177468"
file_view_id <- "syn12179146"
download_id  <- "syn12177447"
index_id     <- "syn12213028"

cell_lines_to_sample <- c("CD4+_T_cells_BGI",
                          "CD8+_T_cells_BGI")

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

# fastq files
manifest_df <- manifest_id %>% 
    create_df_from_synapse_id %>% 
    filter(sample_name %in% cell_lines_to_sample) %>% 
    select(run, sample_name) %>% 
    set_colnames(c("SRR_id", "cell_type"))

file_view_df <- file_view_id %>% 
    str_c("select id, name, parentId from ", .) %>% 
    synTableQuery %>% 
    as.data.frame %>% 
    filter(parentId == download_id) %>% 
    filter(str_detect(name, "fastq.gz$")) %>% 
    mutate(SRR_id = str_sub(name, end = 10)) %>% 
    inner_join(manifest_df) %>% 
    select(id, cell_type)

gz_paths <- llply(file_view_df$id, download_from_synapse, .parallel = T)
commands <- str_c("gunzip ", gz_paths)
l_ply(commands, system, .parallel = T)

file_view_df %>% 
    inset("path", value = unlist(gz_paths)) %>% 
    mutate(path = str_sub(path, end = -4)) %>% 
    write_tsv("fastq_files.tsv")

# kallisto index
index_file <- download_from_synapse(index_id)
new_path <- basename(index_file)
file.rename(index_file, new_path)
system(str_c("gunzip ", new_path))
    
    