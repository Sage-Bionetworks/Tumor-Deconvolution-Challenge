library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(yaml)
library(curl)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/"

cs_id      <- "syn11898291"
mcp_url    <- "http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"
upload_id  <- "syn11898282"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)

activity_obj <- Activity(
    name = "upload",
    description = "create yaml file containing lists of genes used by various deconvolution tools",
    used = list(mcp_url, cs_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/misc/create_deconvolution_gene_lists.R")
)

mcpcounter_genes <- mcp_url %>% 
    curl %>% 
    read_tsv %>% 
    extract2("HUGO symbols") %>% 
    as.list()


synLogin()

cibersort_genes <- cs_id %>% 
    create_df_from_synapse_id %>% 
    extract2("Gene symbol") %>% 
    as.list()


list(
    "mcpcounter_genes" = mcpcounter_genes,
    "cibersort_genes" = cibersort_genes) %>% 
    write_yaml("deonconvolution_genes.yaml")

upload_file_to_synapse("deonconvolution_genes.yaml", upload_id, activity_obj = activity_obj)
