library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn18070253",
    wikiName = "QC", 
    overwrite = F)
