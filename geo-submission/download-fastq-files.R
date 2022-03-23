library(pacman)
library(synapser)
library(plyr)
synLogin()

suppressPackageStartupMessages(p_load(foreach))
suppressPackageStartupMessages(p_load(parallel))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load('doMC'))
  cat(paste0('Registering ', num.cores-1, ' cores.\n'))
  registerDoMC(cores=(num.cores-1))
}


input.folder.synId <- "syn21557721"
children <- synGetChildren(input.folder.synId)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

annotation.keys <- c("pair", "reads", "sample", "size", "type")

synIds <- as.character(df$id)
names(synIds) <- synIds

df$id <- as.character(df$id)

download.dir <- "/fastscratch/whitebr/download/"
dir.create(download.dir, recursive=TRUE)

new.download.dir <- "/fastscratch/whitebr/download-new/"
dir.create(new.download.dir, recursive=TRUE)

l_ply(synIds, .parallel = TRUE,
        .fun = function(synId) {
                 obj <- synGet(synId, downloadFile=FALSE)
                 if(!file.exists(paste0(download.dir, "/", obj$properties$name)) && !file.exists(paste0(new.download.dir, "/", obj$properties$name))) {
                   cat(paste0("Downloading ", obj$properties$name, "\n"))
# NB: downloading to a different download location
                   obj <- synGet(synId, downloadFile=TRUE, downloadLocation=new.download.dir)
                 } else {
                   cat(paste0("File exists ", obj$properties$name, "\n"))
                 }
               })
