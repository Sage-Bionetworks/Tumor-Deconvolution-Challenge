library(synapser)
library(plyr)
synLogin()

input.folder.synId <- "syn21571479"
children <- synGetChildren(input.folder.synId)
l <- as.list(children)
df <- do.call(rbind.data.frame, l)

synIds <- as.character(df$id)
names(synIds) <- synIds

df$id <- as.character(df$id)

stats <-
  ldply(synIds,
        .fun = function(synId) {
print(synId)
                 f <- synGet(synId, downloadFile=FALSE)
                 md5 = f$get("md5")
                 sz = f$get("fileSize")
                 data.frame(md5 = md5, size = sz)
               })
colnames(stats)[1] <- "id"

m <- merge(df[, !(colnames(df) %in% c("type"))], stats)
                 
write.table(m, file = "processed-file-table.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)