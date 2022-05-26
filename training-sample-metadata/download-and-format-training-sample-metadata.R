suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(library(synapser))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(data.table))
synLogin()

download.and.save.xlsx.as.csv <- function(synId) {
  path <- synGet(synId, downloadFile = TRUE, downloadLocation = ".")$path
  tbl <- read.xlsx(path, sheet = 1)

  output.path <- gsub(path, pattern = "xlsx", replacement = "csv")  
  write.table(tbl, file = output.path, row.names = FALSE, col.names = TRUE, quote = TRUE, sep = ",")
  # Ensure we can read the file path without problems. Only fread seems smart enough to figure this out.
#  tbl <- tryCatch(read.table(output.path, sep = ",", header = TRUE, quote = "\""),
#                  error = function(e) stop(paste0("Trouble reading ", output.path, "\n")))
  tbl <- tryCatch(fread(output.path),
                  error = function(e) stop(paste0("Trouble reading ", output.path, "\n")))
}

# Download sample-level GEO array data
synId <- "syn18728081"
download.and.save.xlsx.as.csv(synId)

# Download sample-level GEO rna-seq data
synId <- "syn18751454"
download.and.save.xlsx.as.csv(synId)