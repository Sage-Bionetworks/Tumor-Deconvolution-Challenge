library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(yaml)
library(ggfortify)
library(RColorBrewer)
library(heatmap.plus)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/"

annotation_id <- "syn11898281"
count_id      <- "syn11898217"
yaml_id       <- "syn11898303"

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

gene_yaml <- yaml_id %>% 
    download_from_synapse %>% 
    read_yaml

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

mean_cpm_matrix <- sample_metadata_df %>%
    split(.$combined_name) %>%
    map(use_series, cell_id) %>%
    map(get_summary_by_matrix_cols, cpm_matrix, mean) %>%
    do.call("cbind", .) %>% 
    .[!rowSums(.) == 0,]

mcp_zscore_matrix <- mean_cpm_matrix %>% 
    .[rownames(.) %in% gene_yaml$mcpcounter_genes,] %>% 
    zscore_matrix
    
cs_zscore_matrix <- mean_cpm_matrix %>% 
    .[rownames(.) %in% gene_yaml$cibersort_genes,] %>% 
    zscore_matrix



color_df1 <- data_frame(
    "cell_type" = unique(sample_metadata_df$cell_type),
    "color_cell_type" = brewer.pal(nlevels(as.factor(sample_metadata_df$cell_type)), name = 'Dark2'))

color_df2 <- data_frame(
    "person" = unique(sample_metadata_df$person),
    "color_person" = brewer.pal(nlevels(as.factor(sample_metadata_df$person)), name = 'Set3'))

person_metadata_df <- sample_metadata_df %>% 
    select(-cell_id) %>% 
    distinct %>% 
    arrange(person, cell_type) %>% 
    left_join(color_df1) %>% 
    left_join(color_df2)

color_matrix <- person_metadata_df %>% 
    select(color_cell_type, color_person) %>% 
    set_colnames(c("Cell_Type", "Person")) %>% 
    as.matrix

autoplot(prcomp(t(mcp_zscore_matrix)), data = person_metadata_df, colour = 'person')
autoplot(prcomp(t(mcp_zscore_matrix)), data = person_metadata_df, colour = 'cell_type')
autoplot(prcomp(t(cs_zscore_matrix)), data = person_metadata_df, colour = 'person')
autoplot(prcomp(t(cs_zscore_matrix)), data = person_metadata_df, colour = 'cell_type')




heatmap.plus(mcp_zscore_matrix, ColSideColors = color_matrix, scale = "none")

heatmap.plus(cs_zscore_matrix, ColSideColors = color_matrix, scale = "none")

