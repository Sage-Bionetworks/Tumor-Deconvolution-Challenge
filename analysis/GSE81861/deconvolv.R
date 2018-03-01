library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(MCPcounter)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE81861/"

annotation_id <- "syn11898281"
count_id      <- "syn11898217"

setwd(home_dir)
source("scripts/utils.R")
source("analysis/GSE81861/dataset_utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

annotation_df <- create_annotation_df(annotation_id)

count_df <- create_df_from_synapse_id(count_id) 

sample_metadata_df <- count_df %>% 
    .[1:3] %>% 
    left_join(annotation_df, by = c("cell_id" = "sample")) %>% 
    mutate(cell_type = replace_na(cell_type, "Unknown")) %>% 
    mutate(combined_name = str_c(person, "_", cell_type))
    
cpm_matrix <- count_df %>% 
    select(-c(sample_id, cell_type)) %>% 
    df_to_matrix("cell_id") %>% 
    t %>% 
    apply(2, calculate_cpm)

mean_cpm_matrix <- aggregate_matrix(cpm_matrix, sample_metadata_df, "person", "cell_id")

mcp_result <- MCPcounter.estimate(mean_cpm_matrix, featuresType = "HUGO_symbols")

mean_cpm_matrix  %>% 
    data.frame %>% 
    rownames_to_column("GeneSymbol") %>% 
    write_tsv("./cpm.tsv")

source("/home/aelamb/repos/irwg/iatlas-tool-cibersort/bin/CIBERSORT.R")
cs_result <- CIBERSORT("/home/aelamb/repos/irwg/iatlas-tool-cibersort/sample.references.matrix.txt",
                       "./cpm.tsv", 
                       QN = F)

mcp_result %>% 
    matrix_to_df("cell_type") %>% 
    write_tsv("mcp_results.tsv")

cs_result %>% 
    matrix_to_df("cell_type") %>% 
    write_tsv("cs_results.tsv")

