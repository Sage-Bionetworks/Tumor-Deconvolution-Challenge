library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(pheatmap)
library(preprocessCore)
library(ggfortify)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE72056/"

expr_id        <- "syn12119637"
anno_id        <- "syn12119636"
genes_id       <- "syn11918430"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

annotation_df <- anno_id %>% 
    create_df_from_synapse_id %>% 
    group_by(cell_type) %>% 
    sample_n(15) %>% 
    ungroup

gene_df <-  create_df_from_synapse_id(genes_id)

mcp_genes <- gene_df %>% 
    filter(Method == "mcpcounter") %>%
    use_series("Hugo")

cs_genes <- gene_df %>% 
    filter(Method == "cibersort") %>%
    use_series("Hugo")

expr_matrix <- expr_id %>%
    create_df_from_synapse_id() %>%
    group_by(Hugo) %>% 
    summarise_all(sum) %>% 
    ungroup %>% 
    df_to_matrix("Hugo") %>% 
    .[,colnames(.) %in% annotation_df$sample]

mcp_heatmap_row_df <- gene_df %>% 
    filter(Method == "mcpcounter") %>% 
    filter(Hugo %in% rownames(expr_matrix)) %>% 
    select(-Method) %>% 
    arrange(cell_type) %>% 
    data.frame %>% 
    column_to_rownames("Hugo") %>% 
    set_names("Cell Type")

heatmap_col_df <- annotation_df %>% 
    data.frame %>% 
    column_to_rownames("sample")

expr_matrix <- expr_matrix %>% 
    .[rownames(mcp_heatmap_row_df),] %>% 
    .[,rownames(heatmap_col_df)]


png('GSE72056_mcpcounter_genes_heatmap.png', width = 1200, height = 1200)
pheatmap(
    expr_matrix,
    main = "MCPCounter GSE72056",
    annotation_row = mcp_heatmap_row_df,
    annotation_col = heatmap_col_df,
    cluster_rows = F,
    cluster_cols = F,
    show_colnames = F,
    scale = "none")
dev.off()
