library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn17091410",
    wikiName = "QC", 
    overwrite = F)
