library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn12649844",
    wikiName = "QC", 
    overwrite = F)