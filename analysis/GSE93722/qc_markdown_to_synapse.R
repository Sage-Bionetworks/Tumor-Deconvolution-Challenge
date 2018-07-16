library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn13363373",
    wikiName = "QC", 
    overwrite = F)