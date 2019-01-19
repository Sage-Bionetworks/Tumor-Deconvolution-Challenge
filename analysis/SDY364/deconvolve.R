library(synapser)
library(tidyverse)
library(data.table)
library(magrittr)
library(MCPcounter)

repo_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/SDY364/"
expr_id   <- "syn17974199"
upload_id <- "syn17974650"


source("../../../irwg/iatlas-tool-cibersort/workflow/docker/cibersort/bin/CIBERSORT.R")
lm22 <- "../../../irwg/iatlas-tool-cibersort/LM22.tsv"
source("../../scripts/utils.R")
synLogin()



expr_df <- expr_id %>% 
    create_df_from_synapse_id 

write_tsv(expr_df, "expr.tsv")


CIBERSORT(lm22, "expr.tsv") %>% 
    as.data.frame %>% 
    tibble::rownames_to_column("sample") %>% 
    readr::write_tsv("cibersort_results.tsv")

expr_df %>% 
    df_to_matrix("Hugo") %>% 
    MCPcounter.estimate(featuresType = "HUGO_symbols") %>% 
    write.table("mcpcounter_results.tsv", sep = "\t")

activity_obj <- Activity(
    name = "create",
    description = "create and upload deconvolution results using cibersort and mcpcounter",
    used = list(expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY364/fastq_mixing_CD4_CD8/deconvolve.R")
)

upload_file_to_synapse("cibersort_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("mcpcounter_results.tsv", upload_id, activity_obj = activity_obj)
