library(synapser)
library(tidyverse)
library(magrittr)

tar_id    <- "syn17079944"
upload_id <- "syn17079942"

source("../../scripts/utils.R")
synLogin()

tar_file <- download_from_synapse(tar_id) 

file.copy(tar_file, ".")

system(str_c("tar -xvf ", basename(tar_file)))

gt_df <- "immunoStates_paper_code_v2/deconvolution_ellison_all_methods_all_matrices.RDS" %>% 
    readRDS() 

gt_df %>% 
    filter(Year == "Ellison_Year2011") %>% 
    ggplot(aes(x = Prop, y = Count, color = Cell)) +
    geom_point() +
    facet_grid(rows = vars(alg), cols = vars(mat)) 

gt_df %>% 
    select(-c("Prop", "alg", "mat")) %>% 
    distinct %>% 
    separate(Sample, into = c("Sample", "Sample_year"), sep = " ") %>% 
    write_tsv("Ellision_ground_truth.tsv")

activity_obj <- Activity(
    name = "create",
    description = "process data from tar into tsv",
    used = list(tar_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/Immunostates/get_files_from_tar.R")
)

upload_file_to_synapse("Ellision_ground_truth.tsv", upload_id, activity_obj = activity_obj)