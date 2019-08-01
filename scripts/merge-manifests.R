library(pacman)

suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(synapserutils))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(magrittr))
suppressPackageStartupMessages(p_load(xlsx))
suppressPackageStartupMessages(p_load(reshape2))

source("utils.R")
synLogin()

top.level.folder <- "syn20188706"

children <- synGetChildren(top.level.folder)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)
df <- subset(df, type == "org.sagebionetworks.repo.model.Folder")

datasets <- as.character(df$name)
names(datasets) <- datasets
folder.synIds <- as.list(df$id)
names(folder.synIds) <- datasets
manifests <-
  llply(datasets,
        .fun = function(ds) {
	         synId <- as.character(folder.synIds[[ds]])
                 children <- synGetChildren(synId)
                 l <- as.list(children)
                 df <- do.call(rbind.data.frame, l)
	         file <- paste0(ds, "-metadata.tsv")
		 flag <- df$name == file
		 if(length(which(flag)) != 1) {
		   return(NULL)
		   stop(paste0("Could not find ", file, " amongst: ",
		               paste(df$name, collapse = ", ")))
		 }
		 synId <- as.character(df[flag, "id"])
		 obj <- synGet(synId, downloadFile = TRUE)
		 file <- get.synapse.file.location(obj)
		 manifest <- read.table(file, sep = "\t", header = TRUE,
		                        as.is = TRUE)
                 vec <- manifest$value
		 names(vec) <- manifest$key
                 vec
	       })

exclude <- c("SDY180", "SDY305", "SDY984", "MCP-counter-RNA-seq")
manifests <- manifests[!(names(manifests) %in% exclude)]

cols <- names(manifests[[1]])

## llply(manifests, .fun = function(df) print(length(df)))

combined <- ldply(manifests, .fun = function(vec) vec[cols])

write.table(combined, file = "leaderboard-dataset-metadata.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE,
	    quote = FALSE)