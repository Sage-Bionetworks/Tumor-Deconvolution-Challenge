library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(magrittr)

home_dir      <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir       <- "/home/aelamb/tmp/"
download_path <- "/tmp/mozilla_aelamb0/sanitized_annotations/"
upload_id     <- "syn11867744"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)

registerDoMC(cores = detectCores())

files <- download_path %>% 
    list.files(full.name = T)

activity_obj <- Activity(
    name = "upload",
    description = "upload unziped sanitized annotations after downloading from dropbox",
    used = list("https://www.dropbox.com/s/ky8t1yc17s40qz8/Archive.zip?dl=0"),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE81861/upload_annotation_files_to_synapse.R")
)

synLogin()
l_ply(files, upload_file_to_synapse, upload_id, activity_obj = activity_obj, .parallel = T)



    


