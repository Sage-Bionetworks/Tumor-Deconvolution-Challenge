library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(yaml)
library(curl)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/"

upload_id     <- "syn11867744"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

mpc_df <- data_frame(
    "mcpcounter_cell_lines" = c(
        "B lineage", "CD8 T cells", "Cytotoxic lymphocytes", "Endothelial cells", 
        "Fibroblasts", "Monocytic lineage", "Myeloid dendritic cells", 
        "Neutrophils", "NK cells", "T cells", NA, NA),
    "dataset_cell_lines" = c(
        "Bcell", NA, NA, "Endothelial", "Fibroblast", "Macrophage", NA, NA, 
        NA, "Tcell", "MastCell", "Epithelial")) 

    
cs_df <- data_frame(
    "cibersort_cell_lines" = c(
        "B cells memory", "B cells naive", "Dendritic cells activated",
        "Dendritic cells resting", "Eosinophils", "Macrophages M0", 
        "Macrophages M1", "Macrophages M2", "Mast cells activated", 
        "Mast cells resting", "Monocytes", "Neutrophils", "NK cells activated",
        "NK cells resting", "Plasma cells", "T cells CD4 memory activated",
        "T cells CD4 memory resting", "T cells CD4 naive", "T cells CD8", 
        "T cells follicular helper", "T cells gamma delta", 
        "T cells regulatory (Tregs)", NA, NA, NA),
    "dataset_cell_lines" = c(
        "Bcell", "Bcell", NA, NA, NA, "Macrophage", "Macrophage", "Macrophage",
        "MastCell", "MastCell", "Monocytic lineage", NA, NA, NA, NA, "Tcell",
        "Tcell", "Tcell", "Tcell", "Tcell", "Tcell", "Tcell", "Endothelial", 
        "Fibroblast", "Epithelial"))

write_tsv(mpc_df, "mcpcounter_translation_table.tsv")
write_tsv(cs_df, "cibersort_translation_table.tsv")
upload_file_to_synapse("mcpcounter_translation_table.tsv", upload_id)
upload_file_to_synapse("cibersort_translation_table.tsv", upload_id)
    