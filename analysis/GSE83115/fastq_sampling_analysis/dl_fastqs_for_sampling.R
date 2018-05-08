library(plyr)
library(doMC)
library(tidyverse)
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

cell_lines_to_sample <- c("Monocytes_BGI", 
                          "CD4+_T_cells_BGI",
                          "CD8+_T_cells_BGI")

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

SRR_ids<- manifest_id %>% 
    create_df_from_synapse_id %>% 
    filter(sample_name %in% cell_lines_to_sample) %>% 
    use_series(run)

file_view_id %>% 
    str_c("select id, name, parentId from ", .) %>% 
    synTableQuery %>% 
    as.data.frame %>% 
    filter(parentId == download_id) %>% 
    filter(str_detect(name, "fastq.gz$")) %>% 
    mutate(SRR_id = str_sub(name, end = 10)) %>% 
    filter(SRR_id %in% SRR_ids) %>% 
    use_series(id) %>% 
    llply(download_from_synapse, .parallel = T) %>% 
    str_c("gunzip ", .) %>%
    l_ply(system, .parallel = T)

