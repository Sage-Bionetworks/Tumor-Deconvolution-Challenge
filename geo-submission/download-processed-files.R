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


download.dir <- "/fastscratch/whitebr/download"
dir.create(download.dir, recursive=TRUE)


input.folder.synId <- "syn21571479"
children <- synGetChildren(input.folder.synId)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

synIds <- as.character(df$id)
names(synIds) <- synIds

df$id <- as.character(df$id)

        
l_ply(synIds,
        .fun = function(synId) {
print(synId)
                 synGet(synId, downloadFile=TRUE, downloadLocation = download.dir)
               })
