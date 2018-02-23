library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/"
download_path <- "/tmp/mozilla_aelamb0/sanitized_annotations"
upload_id     <- "syn11867744"

source_url    <-  "https://www.dropbox.com/s/ky8t1yc17s40qz8/Archive.zip?dl=0"



setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())