library(synapser)
library(knit2synapse)
synLogin()

## wiki_synid is defined in setup.R
source("setup.R")
knitfile2synapse(
    "./qc_markdown_flow.Rmd",
    owner = wiki_synid,
    wikiName = "QC", 
    overwrite = F)
