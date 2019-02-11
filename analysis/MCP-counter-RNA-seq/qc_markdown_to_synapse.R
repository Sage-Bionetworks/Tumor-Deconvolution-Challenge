library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn18103327",
    wikiName = "QC", 
    overwrite = F)
