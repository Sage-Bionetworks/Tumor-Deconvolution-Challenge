library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/E-MTAB-2319/"

annotation_url  <- "https://raw.githubusercontent.com/mdozmorov/63_immune_cells/master/data/E-MTAB-2319.sdrf.txt"
cs_id           <- "syn11898291"
mcp_url         <- "http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"
script_url      <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/E-MTAB-2319/create_annotation_tables.R"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()



annotation_df <- fread(annotation_url)
dataset_types <- annotation_df %>% 
    extract2("Factor Value[cell type]") %>% 
    unique %>% 
    sort

cibersort_types <- cs_id %>% 
    create_df_from_synapse_id %>%
    colnames %>% 
    .[-1] %>% 
    unique %>% 
    sort

mcpcounter_types <- mcp_url %>% 
    fread %>% 
    extract2("Cell population") %>% 
    unique %>% 
    sort

mcp_df <- data_frame(
    "mcpcounter_cell_type" = c(
        "B lineage",  "B lineage",  "B lineage", "T cells", "T cells", 
        "T cells", "T cells", "T cells", "T cells", "T cells", "CD8 T cells",
        "T cells", "CD8 T cells", "T cells", "CD8 T cells", "T cells", 
        "Cytotoxic lymphocytes", "Endothelial cells", "Fibroblasts",
        "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils",
        "NK cells"),
    "translated_mpcounter_cell_type" = c(
        "B CD5", "B Memory", "B Naive", "CD4 Central Memory", 
        "CD4 Effector Memory", "CD4 Naive", "CD4 Th1", "CD4 Th17", "CD4 Th2", 
        "CD4 Treg", "CD8 Central Memory", "CD8 Central Memory", 
        "CD8 Effector Memory", "CD8 Effector Memory", "CD8 Naive", "CD8 Naive", 
        rep(NA, 7)))

cibersort_df <- data_frame(
    "cibersort_cell_type" = c(
        NA, 
        "B cells memory", 
        "B cells naive", 
        "T cells CD4 memory activated",
        "T cells CD4 memory activated",
        "T cells CD4 memory resting",
        "T cells CD4 memory resting",
        "T cells CD4 naive",
        "T cells regulatory (Tregs)",
        NA,
        NA,
        NA,
        "T cells follicular helper", 
        "T cells gamma delta",
        "T cells CD8",
        "T cells CD8",
        "T cells CD8",
        "Dendritic cells resting", "Eosinophils", "Macrophages M0",
        "Macrophages M1", "Macrophages M2", "Mast cells activated",
        "Mast cells resting", "Monocytes", "Neutrophils", "NK cells activated",
        "NK cells resting", "Plasma cells", "Dendritic cells activated"
        ),
    "translated_cibersort_cell_type" = c(
        "B CD5", 
        "B Memory", 
        "B Naive", 
        "CD4 Central Memory", 
        "CD4 Effector Memory",
        "CD4 Central Memory", 
        "CD4 Effector Memory",
        "CD4 Naive", 
        "CD4 Treg",
        "CD4 Th1", 
        "CD4 Th17", 
        "CD4 Th2", 
        NA,
        NA,
        "CD8 Central Memory",
        "CD8 Effector Memory", 
        "CD8 Naive", 
        rep(NA, 13)))

write_tsv(annotation_df, "annotation_df.tsv")
write_tsv(mpc_df, "mcpcounter_translation_table.tsv")
write_tsv(cibersort_df, "cibersort_translation_table.tsv")


activity_obj <- Activity(
    name = "upload",
    description = "format annotation table from url",
    used = list(annotation_url),
    executed = list(script_url)
)

upload_file_to_synapse("annotation_df.tsv", "syn11958707", activity_obj = activity_obj)

activity_obj <- Activity(
    name = "upload",
    description = "create file",
    used = list(annotation_url, mcp_url),
    executed = list(script_url)
)

upload_file_to_synapse("mcpcounter_translation_table.tsv", "syn11958705", activity_obj = activity_obj)

activity_obj <- Activity(
    name = "upload",
    description = "create file",
    used = list(annotation_url, cs_id),
    executed = list(script_url)
)

upload_file_to_synapse("cibersort_translation_table.tsv", "syn11958705", activity_obj = activity_obj)