library(synapser)
library(knit2synapse)
synLogin()
knitfile2synapse(
    "./qc_markdown.Rmd",
    owner = "syn17971750",
    wikiName = "QC", 
    overwrite = F)
