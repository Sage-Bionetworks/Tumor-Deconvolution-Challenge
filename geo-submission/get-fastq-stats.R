library(synapser)
library(plyr)
synLogin()

input.folder.synId <- "syn21557721"
children <- synGetChildren(input.folder.synId)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

annotation.keys <- c("pair", "reads", "sample", "size", "type")

synIds <- as.character(df$id)
names(synIds) <- synIds

df$id <- as.character(df$id)

stats <-
  ldply(synIds,
        .fun = function(synId) {
print(synId)
                 f <- synGet(synId, downloadFile=FALSE)
                 md5 = f$get("md5")
                 pair = f$annotations[["pair"]][0]
                 if(is.null(pair)) { pair <- NA }
                 reads = f$annotations[["reads"]][0]
                 if(is.null(reads)) { reads <- NA }
                 smpl = f$annotations[["sample"]][0]
                 if(is.null(smpl)) { smpl <- NA }
                 sz = f$annotations[["size"]][0]
                 if(is.null(sz)) { sz <- NA }
                 typ = f$annotations[["type"]][0]
                 if(is.null(typ)) { typ <- NA }
                 data.frame(md5 = md5, pair = pair, reads = reads, sample = smpl, size = sz, type = typ)
               })
colnames(stats)[1] <- "id"

m <- merge(df[, !(colnames(df) %in% c("type"))], stats)
                 
write.table(m, file = "fastq-table.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
