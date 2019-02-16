library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)

repo_dir  <- "../../../Tumor-Deconvolution-Challenge/"
tmp_dir   <- tempdir()
tmp_dir <- "."

## Define deconvolutoin_results_folder_upload_id and dataset in setup.R
source("setup.R")

## linear_expr_id is GSEXXX/pre-processed/
## This expression file should be in linear/non-log space!
linear_expr_id   <- "syn18345213"
postfix <- "_labs"

github.path <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
dataset.path <- paste0(github.path, dataset)
deconvolve.url <- paste0(dataset.path, "/deconvolve.R")

source(str_c(repo_dir, "scripts/synapseClient_functions.R"))
source(str_c(repo_dir, "scripts/matrix_functions.R"))
source(str_c(repo_dir, "scripts/deconvolution-utils.R"))
setwd(tmp_dir)
synapseLogin()

paths <- download.deconvolution.tools.and.matrices()
cs.cwltool <- paths$cs.cwltool
cs.ref.matrix <- paths$cs.ref.matrix
mcp.cwltool <- paths$mcp.cwltool

run.deconvolution.tools(linear_expr_id, deconvolution_results_folder_upload_id, cs.cwltool, cs.ref.matrix,
                        mcp.cwltool, exec.url = deconvolve.url,
                        cibersort.output.file = paste0("cibersort_results", postfix, ".tsv"),
                        cibersort.abs.output.file = paste0("cibersort_abs_results", postfix, ".tsv"),
                        mcp.counter.output.file = paste0("mcpcounter_results", postfix, ".tsv"))
