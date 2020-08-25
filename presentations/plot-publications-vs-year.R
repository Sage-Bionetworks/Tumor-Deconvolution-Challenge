library(synapser)
library(plyr)
library(ggplot2)

synLogin()

pub.tbl.synId <- "syn18818793"

tsv <- as.data.frame(synTableQuery(paste0("select * from ", pub.tbl.synId)))

tsv <- tsv[order(tsv$Year, decreasing=FALSE),]
tsv$Publication.Year <- as.factor(tsv$Year)

tsv <-
    ddply(tsv, .variables = "Publication.Year",
          .fun = function(df) data.frame(cnt = nrow(df)))
tsv$cum.cnt <- cumsum(tsv$cnt)

g <- ggplot(data = tsv)
g <- g + geom_col(aes(x = Publication.Year, y = cum.cnt))
g <- g + xlab("Publication Year") + ylab("Cumulative Publication Count")
g <- g + theme(axis.text.x = element_text(angle = 90))
g <- g + theme(text = element_text(size = 20))
png("deconvolution-pubs.png")
print(g)
d <- dev.off()

