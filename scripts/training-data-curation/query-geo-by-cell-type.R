#!/usr/bin/env Rscript

library(GEOmetadb)

cell.type.file <- "phenotypes-to-query.tsv"

sfile <- "GEOmetadb.sqlite"

if(!file.exists(sfile)) getSQLiteFile()

library(dplyr)
library(plyr)
library(parallel)

gmdb = src_sqlite(sfile)
## List available tables in the database
##src_tbls(gmdb)

##con <- dbConnect(SQLite(),'GEOmetadb.sqlite')

##dbDisconnect(con)

tgsm = tbl(gmdb, 'gsm')
tgse = tbl(gmdb, 'gse')
tgse.gsm = tbl(gmdb, 'gse_gsm')
tgpl = tbl(gmdb, 'gpl')

tgpl.df <- as.data.frame(tgpl %>% filter(technology == "high-throughput sequencing") %>% collect())
## Limit to illumina sequencing platforms
tgpl.df <- tgpl.df[grepl(pattern="illumina", ignore.case=TRUE, tgpl.df$title),]
tgpl.df <- tgpl.df[grepl(pattern="sapiens", ignore.case=TRUE, tgpl.df$organism),]
## Limit to HiSeq and NextSeq
tgpl.df <- tgpl.df[grepl(pattern="hiseq", ignore.case=TRUE, tgpl.df$title) | grepl(pattern="nextseq", ignore.case=TRUE, tgpl.df$title),]


## tgsm.sub <- subset(tgsm.sub, gpl %in% tgpl.df$gpl)
## tgsm.sub <- merge(tgsm.sub, tgpl.df[, c("gpl", "title")])


## Get all of the sequencing GSEs
technology <- "Expression profiling by high throughput sequencing"
tmp <- tgse %>% filter(type %in% technology) %>% collect()

tgse.pub <- as.data.frame(tmp)[, c("gse", "pubmed_id", "web_link")]

tgse.gsm.sub <- tgse.gsm %>% filter(gse %in% tmp$gse) %>% collect() 

tgsm.sub <- tgsm %>% filter( (gsm %in% tgse.gsm.sub$gsm) & (organism_ch1 == "Homo sapiens") ) %>% collect()
## Limit to illumina sequencing platform
## tgsm.sub <- subset(tgsm.sub, gpl %in% tgpl.df$gpl)
tgpl.df <- tgpl.df[, c("gpl", "title")]
colnames(tgpl.df) <- c("gpl", "platform")
tgsm.sub <- merge(tgsm.sub, tgpl.df)
tgsm.sub$platform <- unlist(lapply(tgsm.sub$platform, function(str) gsub(str, pattern="([^\\(]+)\\(.+", replacement="\\1")))


cell.type.tbl <- read.table(cell.type.file, sep="\t", header=TRUE)
cols <- colnames(tgsm.sub)[grepl(pattern="characteristic", colnames(tgsm.sub)) | grepl(pattern="source_name", colnames(tgsm.sub))]

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

save.image(".Rdata")

nrows <- nrow(cell.type.tbl)
## nrows <- 2
indices <- 1:nrows
names(indices) <- cell.type.tbl$phenotype[1:nrows]
names(cols) <- cols
bools <- llply(cols, .parallel = TRUE,
               .fun = function(col) {
                        str <- tgsm.sub[[col]]
                        lsts <- llply(indices, .parallel = TRUE,
                              .fun = function(i) {
                                       cell.type <- as.character(cell.type.tbl$phenotype[i])
                                       ignore.case <- cell.type.tbl$ignore.case[i]
                                       cat(paste0("Querying ", col, " for ", cell.type, "\n"))
                                       grepl(pattern=cell.type, x=str, ignore.case=ignore.case)
                              })
##                        Reduce('|', lsts)
                        lsts
               })
## flag <- Reduce('|', bools)

gses <- tgsm.sub[["series_id"]]
gsms <- tgsm.sub[["gsm"]]
submissions <- tgsm.sub[["submission_date"]]
platforms <- tgsm.sub[["platform"]]

df <- llply(cols, .parallel = FALSE,
               .fun = function(col) {
                   str <- tgsm.sub[[col]]
                   str <- gsub(str, pattern="\t", replacement=" ")
                        dfs <- llply(indices, .parallel = FALSE,
                              .fun = function(i) {
                                       cell.type <- as.character(cell.type.tbl$phenotype[i])
                                       flag <- bools[[col]][[cell.type]]
if(!any(flag)) { return(NULL) }
                                       sub <- str[flag]
                                       data.frame(cell.type = cell.type, entry = sub, gse = gses[flag], gsm = gsms[flag], submission = submissions[flag], platform = platforms[flag])
                              })
##                        Reduce('|', lsts)
print(head(dfs))
                        Reduce('rbind', dfs)
               })

## remove line (cell line)
## CD34+
## PC
## germinal
## lc.batch
## osteosarcoma


for(nm in names(df)) {
    if(is.null(df[[nm]])) { next }
    df[[nm]] <- merge(df[[nm]], tgse.pub)
}

all.entries <- Reduce('rbind', df)
un <- unique(all.entries[,c("cell.type", "entry")])
un <- un[!duplicated(un$entry),]
un$entry <- gsub(x=un$entry, pattern="\t", replacement=" ")
write.table(file="entries.xls", un, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

anno.entries <- read.table("annotated-entries.tsv", sep="\t", header=TRUE, comment.char="", quote="'")

anno.entries <- anno.entries[anno.entries$keep == "y",]

uniq.pub.list <- list()
for(nm in names(df)) {
    if(is.null(df[[nm]])) {
        uniq.pub.list[[nm]] <- NULL
    } else {
        uniq.pub.list[[nm]] <- unique(df[[nm]][, c("gse", "pubmed_id", "web_link", "submission", "platform")])
    }
}

anno.uniq.pub.list <- list()
for(nm in names(df)) {
    if(is.null(df[[nm]])) {
        anno.uniq.pub.list[[nm]] <- NULL
    } else {
        foo <- df[[nm]][df[[nm]]$entry %in% anno.entries$entry,]
        anno.uniq.pub.list[[nm]] <- unique(foo[, c("gse", "pubmed_id", "web_link", "submission", "platform")])
    }
}

output.pubs <- function(uniq.pub.list, suffix) {

    all.uniq.pubs <- unique(Reduce('rbind', uniq.pub.list))
    all.uniq.pubs$submission.yr <- unlist(lapply(all.uniq.pubs$submission, function(str) gsub(x=str, pattern="([^-]+)-([^-]+)-([^-]+)", replacement="\\1")))
    all.uniq.pubs$submission.mo <- unlist(lapply(all.uniq.pubs$submission, function(str) gsub(x=str, pattern="([^-]+)-([^-]+)-([^-]+)", replacement="\\2")))
    all.uniq.pubs$submission.day <- unlist(lapply(all.uniq.pubs$submission, function(str) gsub(x=str, pattern="([^-]+)-([^-]+)-([^-]+)", replacement="\\3")))
    
    flag <- !is.na(all.uniq.pubs$pubmed_id)
##    write.table(file="geo-immune-sequencing-pubmed-ids.xls", all.uniq.pubs[flag,],
##                sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    non.nas <- all.uniq.pubs[flag,]
    nas <- all.uniq.pubs[!flag,]
    non.nas <- non.nas[order(non.nas$submission.yr, non.nas$submission.mo, non.nas$submission.day, decreasing=TRUE),]
    foo <- rbind(non.nas, nas)
    write.table(file=paste0("geo-immune-sequencing-pubmed-ids-all", suffix, ".xls"), foo,
                sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

output.pubs(anno.uniq.pub.list, "-anno")

save.image(".Rdata")
