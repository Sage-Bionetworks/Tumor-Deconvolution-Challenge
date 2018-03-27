library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/E-MTAB-2319/"

annotation_url  <- "https://raw.githubusercontent.com/mdozmorov/63_immune_cells/master/data/E-MTAB-2319.sdrf.txt"
script_url      <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/E-MTAB-2319/create_processed_tables.R"
count_id        <- "syn11958709"
hugo_id         <- "syn11536071"
upload_id       <- "syn11958707"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

hugo_df <-  create_df_from_synapse_id(hugo_id)

annotation_df <- annotation_url %>% 
    fread(select = c("Comment[ENA_RUN]", "Factor Value[cell type]")) %>% 
    as_data_frame %>% 
    set_colnames(c("sample", "cell_type")) %>% 
    mutate(sample = str_c(sample, ".bam")) %>% 
    distinct %>% 
    arrange(sample)

counts_df <- count_id %>%
    create_df_from_synapse_id(unzip = T) %>% 
    select(-c(Chr, Start, End, Strand, Length)) %>% 
    dplyr::rename("ensembl_gene_id" = Geneid) %>% 
    left_join(hugo_df) %>% 
    .[,order(colnames(.))] %>% 
    select(ensembl_gene_id, hgnc_symbol, everything()) %>% 
    arrange(ensembl_gene_id)



write_tsv(annotation_df, "annotation_df.tsv")
write_tsv(counts_df, "counts_df.tsv")


activity_obj <- Activity(
    name = "process files",
    description = "process raw files",
    used = list(annotation_url, count_id, hugo_id),
    executed = list(script_url)
)

upload_file_to_synapse("annotation_df.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("counts_df.tsv", upload_id, activity_obj = activity_obj)
