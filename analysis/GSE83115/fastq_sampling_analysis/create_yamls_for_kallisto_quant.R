library(tidyverse)
library(synapser)
library(doMC)

# local
# home_dir     <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# kallisto_dir <- "/home/aelamb/repos/kallisto_cwl/"
# work_dir     <- "/home/aelamb//tmp/tumor_deconvolution/GSE83115/"

# ec2
home_dir     <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
kallisto_dir <- "/home/ubuntu/kallisto_cwl/"
work_dir     <- "/home/ubuntu/"

fastq_file <- "fastq.tsv"
index_id   <- "syn12213028"

setwd(home_dir)
source("scripts/utils.R")
source(str_c(kallisto_dir, "utils/write_yaml.R"))
setwd(work_dir)
synLogin()



manifest_df <- read_tsv(fastq_file) 

df1 <- manifest_df %>% 
    select(-c(name, paths)) %>% 
    spread(key = pair, value = id) %>% 
    set_names(c("sample_name", "p1_id", "p2_id"))

df2 <- manifest_df %>% 
    select(-c(name, id)) %>% 
    spread(key = pair, value = paths) %>% 
    set_names(c("sample_name", "p1_path", "p2_path"))

manifest_df <- left_join(df1, df2) %>% 
    mutate("yaml" = str_c(sample_name, ".yaml")) %>% 
    mutate("log" = str_c(sample_name, ".log")) 

n_cores    <- detectCores() - 1
index_file <- download_from_synapse(index_id)

create_yaml_per_row <- function(row){
    create_kallisto_quant_yaml(
        yaml_file = row$yaml, 
        fastq1 = row$p1_path, 
        fastq2 = row$p2_path, 
        index = index_file,
        threads = n_cores)
}

manifest_df %>% 
    split(1:nrow(.)) %>% 
    walk(create_yaml_per_row)

write_tsv(manifest_df, fastq_file)
