library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(MCPcounter)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE62408/"

count_id      <- "syn11915424"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

rpkm_matrix <- create_df_from_synapse_id(count_id, unzip = T) %>% 
    rename("Hugo" = `-`) %>% 
    group_by(Hugo) %>% 
    summarise_all(mean) %>% 
    df_to_matrix("Hugo")

mcp_result <- MCPcounter.estimate(rpkm_matrix, featuresType = "HUGO_symbols")

rpkm_matrix %>% 
    data.frame %>% 
    rownames_to_column("GeneSymbol") %>% 
    write_tsv("./rpkm.tsv")

source("/home/aelamb/repos/irwg/iatlas-tool-cibersort/bin/CIBERSORT.R")
cs_result <- CIBERSORT("/home/aelamb/repos/irwg/iatlas-tool-cibersort/sample.references.matrix.txt",
                       "./rpkm.tsv",
                       QN = F)

mcp_result %>% 
    matrix_to_df("cell_type") %>% 
    write_tsv("mcp_results.tsv")

cs_result %>% 
    matrix_to_df("cell_type") %>% 
    write_tsv("cs_results.tsv")

activity_obj <- Activity(
    name = "create",
    description = "create and upload deconvolution results using cibersort and mcpcounter cwl files",
    used = list(count_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE62408/deconvolve.R")
)

upload_file_to_synapse("cs_results.tsv", "syn11968719", activity_obj = activity_obj)
upload_file_to_synapse("mcp_results.tsv", "syn11968719", activity_obj = activity_obj)

