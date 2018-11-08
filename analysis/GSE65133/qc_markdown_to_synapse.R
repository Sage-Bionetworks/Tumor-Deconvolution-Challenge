library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn15664931",
    wikiName = "QC", 
    overwrite = F)
