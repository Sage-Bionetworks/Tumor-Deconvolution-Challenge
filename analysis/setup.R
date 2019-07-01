usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, repos = "http://cran.us.r-project.org", dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("pacman")

## Install synapser, if not installed.
## On Ubuntu, this required that I first install ffi, which I did via:
## sudo apt-get install build-essential libssl-dev libffi-dev python-dev
if (!is.element("synapser", installed.packages()[,1])) {
  install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
}

if (!is.element("synapserutils", installed.packages()[,1])) {
  install.packages("synapserutils", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
}

if (!is.element("MCPcounter", installed.packages()[,1])) {
  suppressPackageStartupMessages(p_load("devtools"))
  suppressPackageStartupMessages(p_load("curl"))
  install_github("ebecht/MCPcounter",ref="master", subdir="Source")
}

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

url.base <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"

cibersort.script <- "/home/bwhite/tumor-deconvolution-external/CIBERSORT.R"
cibersort.lm22.signature.file <- "/home/bwhite/tumor-deconvolution-external/LM22.txt"

