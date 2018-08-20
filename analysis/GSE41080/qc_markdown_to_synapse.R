library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn14566979",
    wikiName = "QC", 
    overwrite = F)