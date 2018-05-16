library(tidyverse)
library(synapser)
library(doMC)

# local
# home_dir     <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# workflow_dir <- "/home/aelamb/repos/fastq_mixing_workflow_cwl/"
# work_dir     <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"

# ec2
home_dir     <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
workflow_dir <- "/home/ubuntu/fastq_mixing_workflow_cwl/"
work_dir     <- "/home/ubuntu/"

fastq_file <- "fastq_files.tsv"
index_file <- "GRCH38.idx"
replicates <- 3
CD8_fractions <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50)



setwd(home_dir)
source("scripts/utils.R")
source(str_c(workflow_dir, "utils/write_yaml.R"))
setwd(work_dir)
synLogin()
n_cores <- detectCores() - 1

fastq_df <- fastq_file %>%
    read_tsv %>%
    mutate(pair = str_match(path, "_([12]).fastq$")[,2]) %>% 
    mutate(col = str_c(cell_type, "_", pair)) %>% 
    mutate(col = str_remove_all(col, "\\+")) %>% 
    select(-c(cell_type, pair))

path_df <- fastq_df %>%
    select(- id) %>%
    mutate(col = str_c(col, "_path")) %>%
    spread(key = col, value = path)

yaml_df <-
    data_frame("CD8_fractions" = CD8_fractions) %>%
    mutate(CD4_fractions = 1 - CD8_fractions) %>% 
    mutate(prefix = str_c("CD4_CD8_", CD8_fractions)) %>% 
    mutate(yaml = str_c(prefix, ".yaml")) %>% 
    mutate(log = str_c(prefix, ".log")) %>% 
    inset("seed1", value = sample(1:10000, 7)) %>% 
    inset("seed2", value = sample(1:10000, 7)) %>% 
    inset("seed3", value = sample(1:10000, 7)) %>% 
    merge(path_df) 

create_yaml_by_row <- function(row){
    create_repeat_workflow_yaml(
        row$yaml,
        str_c(row$prefix, "_rep", 1:replicates, ".tsv"),
        c(row$CD4_T_cells_BGI_1_path, 
          row$CD8_T_cells_BGI_1_path),
        c(row$CD4_T_cells_BGI_2_path, 
          row$CD8_T_cells_BGI_2_path),
        c(row$CD4_fractions, 
          row$CD8_fractions),
        c(row$seed1,
          row$seed2,
          row$seed3),
        index_file,
        as.integer(n_cores))
}

yaml_df %>% 
    split(1:nrow(.)) %>% 
    walk(create_yaml_by_row)

yaml_df %>% 
    select(CD8_fractions, CD4_fractions, yaml, log, seed1, seed2, seed3) %>% 
    write_tsv("yaml.tsv")

