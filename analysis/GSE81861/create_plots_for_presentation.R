library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(RColorBrewer)
library(ggfortify)
library(pheatmap)
library(preprocessCore)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE81861/"

annotation_id <- "syn11898281"
count_id      <- "syn11898217"
genes_id       <- "syn11918430"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())


annotation_df <- annotation_id %>% 
    create_df_from_synapse_id("./") %>% 
    select(title, characteristics_ch1) %>% 
    separate(characteristics_ch1, c("string", "person"), ": ") %>% 
    mutate(person = str_replace_all(person, ";", "")) %>% 
    select(title, person) %>% 
    set_names(c("sample", "person"))

count_df <- create_df_from_synapse_id(count_id, "./") 

gene_df <-  create_df_from_synapse_id(genes_id)

mcp_genes <- gene_df %>% 
    filter(Method == "mcpcounter") %>%
    use_series("Hugo")

cs_genes <- gene_df %>% 
    filter(Method == "cibersort") %>%
    use_series("Hugo")

sample_metadata_df <- count_df %>% 
    .[1:3] %>% 
    left_join(annotation_df, by = c("cell_id" = "sample")) %>% 
    mutate(cell_type = replace_na(cell_type, "Unknown")) %>% 
    mutate(combined_name = str_c(person, "_", cell_type))


heatmap_matrix <- count_df %>% 
    select(-c(sample_id, cell_type)) %>% 
    df_to_matrix("cell_id") %>% 
    .[,colSums(.) > 0] %>%  
    t %>% 
    add(1) %>% 
    apply(2, calculate_cpm) %>% 
    log10 %>% 
    quantile_normalize_matrix %>% 
    zscore_matrix

mcp_heatmap_row_df <- gene_df %>% 
    filter(Method == "mcpcounter") %>% 
    filter(Hugo %in% mcp_genes) %>% 
    filter(Hugo %in% rownames(heatmap_matrix)) %>% 
    select(-Method) %>% 
    arrange(cell_type) %>% 
    data.frame %>% 
    column_to_rownames("Hugo") %>% 
    set_names("Reference cell type") 
    
mcp_heatmap_matrix <- heatmap_matrix[rownames(mcp_heatmap_row_df),]

cs_heatmap_matrix <- heatmap_matrix[rownames(heatmap_matrix) %in% cs_genes,]
    
heatmap_col_df <- sample_metadata_df %>%
    select(cell_id, person, cell_type) %>% 
    distinct %>% 
    arrange(cell_id) %>% 
    data.frame %>% 
    column_to_rownames("cell_id") %>% 
    set_names(c("Sample", "Mixture cell type"))

png('GSE81861_mcpcounter_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_heatmap_matrix, 
    main = "MCPCounter genes", 
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    cluster_rows = F, 
    scale = "none")
dev.off()

png('GSE81861_mcpcounter_genes_rows_clustered_heatmap.png', width = 4000, height = 4000)
pheatmap(
    mcp_heatmap_matrix, 
    main = "MCPCounter genes", 
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    scale = "none")
dev.off()

png('GSE81861_cibersort_genes_heatmap.png', width = 4000, height = 4000)
pheatmap(
    cs_heatmap_matrix,
    main = "Cibersort genes", 
    annotation_col = heatmap_col_df,
    scale = "none")
dev.off()


# -----------------------------------------------------------------------------


translation_df <- data_frame(
    "cibersort_cell_type" = c(
        "B cells memory", "B cells naive", "Dendritic cells activated",
        "Dendritic cells resting", "Eosinophils", "Macrophages M0", 
        "Macrophages M1", "Macrophages M2", "Mast cells activated", 
        "Mast cells resting", "Monocytes", "Neutrophils", "NK cells activated",
        "NK cells resting", "Plasma cells", "T cells CD4 memory activated",
        "T cells CD4 memory resting", "T cells CD4 naive", "T cells CD8", 
        "T cells follicular helper", "T cells gamma delta", 
        "T cells regulatory (Tregs)", NA, NA, NA),
    "dataset_cell_type" = c(
        "Bcell", "Bcell", NA, NA, NA, "Macrophage", "Macrophage", "Macrophage",
        "MastCell", "MastCell", NA, NA, NA, NA, NA, "Tcell",
        "Tcell", "Tcell", "Tcell", "Tcell", "Tcell", "Tcell", "Endothelial", 
        "Fibroblast", "Epithelial"))

annotated_fraction_df <- sample_metadata_df %>% 
    split(.$person) %>% 
    map(group_by, combined_name) %>% 
    map(dplyr::summarise, count = n()) %>% 
    map(mutate, annotated_fraction = count / sum(count)) %>% 
    bind_rows %>% 
    select(combined_name, annotated_fraction)
    

cs_result_df <- "cs_results.tsv" %>%
    str_c(tmp_dir, .) %>%
    read_tsv %>%
    transpose_df("cell_type", "cell_type") %>% 
    gather("sample", "result_percentage", CRC01:CRC11) %>% 
    mutate(cell_type = str_replace_all(cell_type, "\\.", " ")) %>% 
    left_join(translation_df, by = c("cell_type" = "cibersort_cell_type")) %>% 
    .[complete.cases(.),] %>% 
    split(.$sample) %>% 
    map(group_by, sample, dataset_cell_type) %>% 
    map(dplyr::summarise, predicted_fraction = sum(result_percentage)) %>% 
    bind_rows %>% 
    mutate(combined_name = str_c(sample, dataset_cell_type, sep = "_")) %>% 
    ungroup

plot_df <- cs_result_df %>% 
    inner_join(annotated_fraction_df) 

png('GSE81861_cibersort_scatterplot.png')
ggplot(plot_df) +
    geom_point(aes(x = annotated_fraction, y = predicted_fraction, color = sample, shape = dataset_cell_type)) +
    xlab("Annotated fraction") +
    ylab("Cibersort predicted fraction") +
    theme_bw() 
dev.off()
