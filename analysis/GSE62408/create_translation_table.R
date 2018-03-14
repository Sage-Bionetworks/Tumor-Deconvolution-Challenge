library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/"

cs_id      <- "syn11898291"
mcp_url    <- "http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"

upload_id  <- "syn11958138"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)

activity_obj <- Activity(
    name = "upload",
    description = "create yaml file containing lists of genes used by various deconvolution tools",
    used = list(mcp_url, cs_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE62408/create_translation_table.R")
)

synLogin()

mpc_df <- data_frame(
    "mcpcounter_cell_type" = c(
        "B lineage", "CD8 T cells", "Cytotoxic lymphocytes", "Endothelial cells", 
        "Fibroblasts", "Monocytic lineage", "Myeloid dendritic cells", 
        "Neutrophils", "NK cells", "T cells", "T cells", "T cells"),
    "translated_mpcounter_cell_type" = c(
        "Bcells", "CD8Tcells", NA, NA, NA, "Monocytes", NA, "Granulocytes", 
        "NKcells", "MemoryTcells", "CD4Tcells", "CD8Tcells"))
    
cibersort_df <- data_frame(
    "cibersort_cell_type" = c(
        "B cells memory", "B cells naive", "Dendritic cells activated",
        "Dendritic cells resting", "Eosinophils", "Macrophages M0",
        "Macrophages M1", "Macrophages M2", "Mast cells activated",
        "Mast cells resting", "Monocytes", "Neutrophils", "NK cells activated",
        "NK cells resting", "Plasma cells", "T cells CD4 memory activated",
        "T cells CD4 memory activated", "T cells CD4 memory resting",
        "T cells CD4 memory resting", "T cells CD4 naive", "T cells CD8",
        "T cells follicular helper", "T cells gamma delta",
        "T cells regulatory (Tregs)"),
    "translated_cibersort_cell_type" = c(
        "Bcells", "Bcells", NA, NA, "Granulocytes", NA, NA, NA, "Granulocytes",
        "Granulocytes", "Monocytes", "Granulocytes", "NKcells", "NKcells", NA,
        "MemoryTcells", "CD4Tcells", "MemoryTcells", "CD4Tcells", "CD4Tcells",
        "CD8Tcells", "CD4Tcells", NA, "CD4Tcells"))

write_tsv(mpc_df, "mcpcounter_translation_table")
write_tsv(cibersort_df, "cibersort_translation_table")
upload_file_to_synapse("mcpcounter_translation_table", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("cibersort_translation_table", upload_id, activity_obj = activity_obj)
    