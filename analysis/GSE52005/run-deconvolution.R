source("dataset-setup.R")
source("../deconvolution-methods/mcpcounter-coarse-grained.R")
source("../deconvolution-methods/cibersort-models.R")
source("../deconvolution-methods/deconvolution-models.R")

suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))

## Begin configuration

file <- "run-deconvolution.R"

## End configuration

script_url <- paste0(url.base, "/", dataset, "/", file)

source("../../scripts/utils.R")
synLogin()

output.folder.synId <- get.synapse.id.of.folder.content(staging.folder.synId, dataset, exit.on.failure = TRUE)

metadata.file.name <- paste0(dataset, "-metadata.tsv")

metadata.synId <- get.synapse.id.of.folder.content(output.folder.synId, metadata.file.name, exit.on.failure = TRUE)

run.deconvolution.methods(metadata.synId, output.folder.synId) 

