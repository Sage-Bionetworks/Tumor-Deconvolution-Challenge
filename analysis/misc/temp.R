library(data.table)
library(magrittr)
library(tidyverse)

studies <- c("SDY212", "SDY311", "SDY312", "SDY314", "SDY315", "SDY404")

setwd("/home/aelamb/Downloads/SDY")

cytof_pbmc_df <- "10KImmunomes.CyTOF_ PBMC.2018-06-25.csv" %>% 
    fread %>% 
    as_data_frame %>% 
    filter(study_accession %in% studies)

flow_pbmc_df <- "10KImmunomes.Flow Cytometry_ PBMC.2018-06-25.csv" %>% 
    fread %>% 
    as_data_frame %>% 
    filter(study_accession %in% studies)

flow_whole_blood_df <- "10KImmunomes.Flow Cytometry_ Whole Blood.2018-06-25.csv" %>% 
    fread %>% 
    as_data_frame %>% 
    filter(study_accession %in% studies)

cell_type_df <- bind_rows(cytof_pbmc_df, flow_pbmc_df, flow_whole_blood_df)

cell_type_anno_df <- select(cell_type_df , study_accession, age, gender, race, subject_accession)

cell_type_df <- select(cell_type_df, -c(study_accession, age, gender, race))


expr_pbmc_df <- "10KImmunomes.Gene Expression_ PBMC.2018-06-25.csv" %>% 
    fread %>%
    data.frame %>% 
    filter(study_accession %in% studies)

anno_pbmc_df <- expr_pbmc_df %>% 
    as_data_frame %>% 
    select(subject_accession, data_accession, age, gender, race, study_accession)

expr_pbmc_df <- expr_pbmc_df %>% 
    select(-c(data_accession, age, gender, race, study_accession)) %>% 
    column_to_rownames("subject_accession") %>% 
    as.matrix %>% 
    t %>% 
    as.data.frame %>% 
    rownames_to_column("Hugo") %>% 
    as_data_frame

expr_wb_df <- "10KImmunomes.Gene Expression_ Whole Blood.2018-06-25.csv" %>% 
    fread %>%
    data.frame %>% 
    filter(study_accession %in% studies)

anno_wb_df <- expr_wb_df %>% 
    as_data_frame %>% 
    select(subject_accession, data_accession, age, gender, race, study_accession)

expr_wb_df <- expr_wb_df %>% 
    select(-c(data_accession, age, gender, race, study_accession)) %>% 
    column_to_rownames("subject_accession") %>% 
    as.matrix %>% 
    t %>% 
    as.data.frame %>% 
    rownames_to_column("Hugo") %>% 
    as_data_frame

expr_anno_df <- bind_rows(anno_pbmc_df, anno_wb_df)



