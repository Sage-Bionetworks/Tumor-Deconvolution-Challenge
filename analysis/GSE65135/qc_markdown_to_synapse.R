library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn15664985",
    wikiName = "QC", 
    overwrite = F)
