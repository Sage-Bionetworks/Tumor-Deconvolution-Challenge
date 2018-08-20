library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn15588640",
    wikiName = "QC", 
    overwrite = F)