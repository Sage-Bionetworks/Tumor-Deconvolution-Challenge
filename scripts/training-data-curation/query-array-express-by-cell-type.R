#!/usr/bin/env Rscript

library(ArrayExpress)

## Query ArrayExpress website and filter by sequencing assay and
## homo sapiens.  Save in a file and input to this script.
## Currently in (downloaded on March 23, 2018)
## ArrayExpression-homo-sapiens-sequencing-query-032318.tsv


queryAE <-
function (keywords = NULL, species = NULL, method = "wget") 
{
    if (is.null(keywords) && is.null(species)) 
        stop("No keywords or species specified")
    baseURL = "http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments"
    qr = paste(baseURL, "?keywords=", keywords, "&species=", 
        species, sep = "")
    qr = URLencode(qr)
    queryfilename = paste("query", keywords, species, ".xml", 
        sep = "")
    query = try(download.file(qr, queryfilename, mode = "wb", method = method))
    x = XML::xmlTreeParse(queryfilename)
    ID = sapply(1:length(XML::xmlRoot(x)), function(i) unlist(XML::xmlElementsByTagName(XML::xmlRoot(x)[[i]], 
        "accession"))["accession.children.text.value"])
    names(ID) = NULL
    x2 = XML::xmlTreeParse(queryfilename, useInternalNodes = TRUE)
    ra = XML::getNodeSet(x2, "/experiments//raw[@count]")
    Rawids = sapply(ra, function(r) XML::xmlGetAttr(r, "name"))
    Rawids = gsub(".raw.*.zip", "", Rawids)
    Raw = rep(NA, length(ID))
    names(Raw) = ID
    Raw[which(ID %in% Rawids)] = "yes"
    Raw[which(!ID %in% Rawids)] = "no"
    pr = XML::getNodeSet(x2, "/experiments//fgem[@count]")
    Procids = sapply(pr, function(r) XML::xmlGetAttr(r, "name"))
    Procids = gsub(".processed.*.zip", "", Procids)
    Processed = rep(NA, length(ID))
    names(Processed) = ID
    Processed[which(ID %in% Procids)] = "yes"
    Processed[which(!ID %in% Procids)] = "no"
    date = ArrayExpress:::getelt(x, node = "releasedate", element = "releasedate.children.text.value")
    pmid = ArrayExpress:::getelt(x, node = "bibliography", element = "bibliography.children.accession.children.text.value")
    spec = ArrayExpress:::getelt(x, node = "species", element = "species.children.text.value")
    experimentdesign = ArrayExpress:::getelt(x, node = "experimentdesign", element = "experimentdesign.children.text.value")
    sampleattribute <- ArrayExpress:::geteltmulti(x, node = "sampleattribute",  element1 = "children.category.children.text.value", element2 = "children.value.children.text.value")

    experimentalfactor = ArrayExpress:::geteltmulti(x, node = "experimentalfactor", 
        element1 = "children.name.children.text.value", element2 = "children.value.children.text.value")
    xmlparsed = data.frame(ID = ID, Raw = Raw[ID], Processed = Processed[ID], 
        ReleaseDate = date, PubmedID = pmid, Species = spec, 
        ExperimentDesign = experimentdesign, ExperimentFactors = experimentalfactor)
    return(xmlparsed)
}

## sets = queryAE(keywords = "sequencing assay", species = "homo+sapiens")
## eset = queryAE(keywords = "E-MTAB-2319", species = "homo+sapiens")

library(XML)

parse.array.express.xml <- function(xml.file) {
    data <- xmlParse(xml.file)
    xml_data <- xmlToList(data)

    ## Entries/columns to keep
    cols <- c("species", "sampleattribute", "experimentalfactor", "accession")
    xml_data <- xml_data$experiment[names(xml_data$experiment) %in% cols]

    ## sampleattribute's have a single category with value's
    ## experimentalfactor's have a single name with value's
    for(i in 1:length(xml_data)) {
        if(names(xml_data)[i] == "sampleattribute") {
            nm <- unname(xml_data[i]$sampleattribute[1])
            vals <- unlist(unname(xml_data[i]$sampleattribute[2:length(xml_data[i]$sampleattribute)]))
            xml_data[[i]] <- vals
            names(xml_data)[i] <- nm
        }
        if(names(xml_data)[i] == "experimentalfactor") {
            nm <- unname(xml_data[i]$experimentalfactor[1])
            vals <- unlist(unname(xml_data[i]$experimentalfactor[2:length(xml_data[i]$experimentalfactor)]))
            xml_data[[i]] <- vals
            names(xml_data)[i] <- nm
        }
    }
    as.data.frame(xml_data)
}

query.array.express.keyword <-
function (keywords = NULL, species = NULL, method = "wget") 
{
    if (is.null(keywords) && is.null(species)) 
        stop("No keywords or species specified")
    baseURL = "http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments"
    qr = paste(baseURL, "?keywords=", keywords, "&species=", 
        species, sep = "")
    qr = URLencode(qr)
    queryfilename = paste("query", keywords, species, ".xml", 
        sep = "")
    query = try(download.file(qr, queryfilename, mode = "wb", method = method))
    parse.array.express.xml(queryfilename)
}

meta <- read.table("ArrayExpression-homo-sapiens-sequencing-query-032318.tsv", sep="\t", header=TRUE, as.is = TRUE, comment.char="", quote="\"")

indices <- 1:nrow(meta)
names(indices) <- meta$Accession

## names(accession.ids) <- accession.ids
## ae.dfs <- llply(accession.ids[1:5], .fun = function(id) query.array.express.keyword(keywords = id))

ae.dfs <- llply(indices, .parallel = FALSE,
                .fun = function(indx) {
                    id <- meta$Accession[indx]
                    cat(paste0("\n\nProcessing indx = ", indx, " id = ", id, "\n"))
                    url <- meta$ArrayExpress.URL[indx]
                    ofile <- paste0(id, ".sdrf.txt")
                    if(!file.exists(ofile) || (file.size(ofile) == 0)) {
                        file <- paste0(gsub(x = url, pattern = "experiments", replacement = "files"), ofile)
                        method <- "wget"
                        query = tryCatch(download.file(file, ofile, mode = "wb", method = method),
                                         error = function(e) { return(NULL) })
                        if(is.null(query)) { return(NULL) }
                    }
                    df <- tryCatch(read.table(ofile, sep="\t", header=TRUE),
                                   error = function(e) { return(NULL) })
                })

save.image(".Rdata.ae.dfs")

ae.idfs <- llply(indices, .parallel = FALSE,
                .fun = function(indx) {
                    id <- meta$Accession[indx]
                    cat(paste0("\n\nProcessing indx = ", indx, " id = ", id, "\n"))
                    url <- meta$ArrayExpress.URL[indx]
                    ofile <- paste0(id, ".idf.txt")
                    if(!file.exists(ofile) || (file.size(ofile) == 0)) {
                        file <- paste0(gsub(x = url, pattern = "experiments", replacement = "files"), ofile)
                        method <- "wget"
                        query = tryCatch(download.file(file, ofile, mode = "wb", method = method),
                                         error = function(e) { return(NULL) })
                        if(is.null(query)) { return(NULL) }
                    }
                    df <- tryCatch(readLines(ofile),
                                   error = function(e) { return(NULL) })
                })

save.image(".Rdata.ae.idfs")

## Ensure that our E-MTAB-2319 set of interest is included here
any(names(ae.dfs) == "E-MTAB-2319")

## sets = queryAE(keywords = "pneumonia", species = "homo+sapiens")

library(dplyr)
library(plyr)
library(parallel)

cell.type.file <- "phenotypes-to-query.tsv"
cell.type.tbl <- read.table(cell.type.file, sep="\t", header=TRUE)

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}


nrows <- nrow(cell.type.tbl)
indices <- 1:nrows
names(indices) <- cell.type.tbl$phenotype[1:nrows]
df.indices <- 1:length(ae.dfs)
names(df.indices) <- names(ae.dfs)
## df.indices <- df.indices[c(1,3930)]
processed.dfs <- 
    llply(df.indices, .parallel = TRUE,
          .fun = function(df.indx) {
              nm <- names(ae.dfs)[df.indx]
              df <- ae.dfs[[nm]]
              cat(paste0("Processing df ", df.indx, ": ", nm, "\n"))
              cols <- colnames(df)
              cols <- cols[cols != "Source.Name"]
              names(cols) <- cols
              dfs2 <- llply(cols, .parallel = FALSE,
                             .fun = function(col) {
                                 str <- df[,col]
                                 dfs <- llply(indices, .parallel = TRUE,
                                               .fun = function(i) {
                                                   cell.type <- as.character(cell.type.tbl$phenotype[i])
                                                   ignore.case <- cell.type.tbl$ignore.case[i]
                                                   ## cat(paste0("Querying ", col, " for ", cell.type, "\n"))
                                                   flag <- grepl(pattern=cell.type, x=str, ignore.case=ignore.case)
                                                   if(!any(flag)) { return(NULL) }
                                                   sub <- str[flag]
                                                   data.frame(cell.type = cell.type, entry = sub, col = col, Source.Name = df$Source.Name[flag])
                                               })
                                 Reduce('rbind', dfs)
                             })
              Reduce('rbind', dfs2)              
          })

save.image(".Rdata.processed.dfs")

null.flag <- unlist(lapply(processed.dfs, is.null))
non.null.processed.dfs <- processed.dfs[!null.flag]

full.df <- ldply(non.null.processed.dfs, data.frame)
colnames(full.df)[1] <- "Accession"

full.df <- unique(full.df)

condensed.df <- ddply(full.df, .variables = c("Accession", "cell.type", "Source.Name"),
                      .fun = function(df) {
                          entry <- paste(df$entry, collapse=",")
                          col <- paste(df$col, collapse=",")
                          data.frame(entry = entry, col = col)
                      })

num.samples.by.accession.cell.type.col <-
    ddply(condensed.df, .variables = c("Accession", "cell.type", "entry", "col"),
          .fun = function(df) {
              data.frame(num.samples = nrow(df))
          })

## Unique by accession and cell type
cell.type.by.accession <- unique(condensed.df[, c("Accession", "cell.type", "entry", "col")])

## Let's exclude PC (plasma cells), DC (dendritic cells), and GC (germincal center)
cell.type.by.accession <- cell.type.by.accession[!(cell.type.by.accession$cell.type %in% c("PC", "DC", "GC")),]
cell.type.by.accession$cell.type <- as.character(cell.type.by.accession$cell.type)

## Merge Accession and entry -- only keep one per entry

num.cell.types.by.accession <- 
    ddply(unique(cell.type.by.accession[, c("Accession", "cell.type")]), .variables = c("Accession"),
          .fun = function(df) nrow(df))
colnames(num.cell.types.by.accession)[2] <- "num.cell.types"

cell.type.by.accession <- merge(cell.type.by.accession, num.cell.types.by.accession)

cell.type.by.accession <- cell.type.by.accession[order(cell.type.by.accession$num.cell.types, cell.type.by.accession$Accession, decreasing = TRUE),]

flag <- !duplicated(cell.type.by.accession$Accession)
cell.type.by.accession$keep <- rep("", nrow(cell.type.by.accession))
cell.type.by.accession$keep[flag] <- FALSE

cell.type.by.accession.with.sample.counts <- merge(cell.type.by.accession, num.samples.by.accession.cell.type.col)

cell.type.by.accession.with.sample.counts <- cell.type.by.accession.with.sample.counts[order(cell.type.by.accession.with.sample.counts$num.cell.types, cell.type.by.accession.with.sample.counts$Accession, cell.type.by.accession.with.sample.counts$num.samples, decreasing = TRUE),]

flag <- !duplicated(cell.type.by.accession.with.sample.counts$Accession)
cell.type.by.accession.with.sample.counts$keep <- rep("", nrow(cell.type.by.accession.with.sample.counts))
cell.type.by.accession.with.sample.counts$keep[flag] <- FALSE


write.table(file = "cell-type-by-accession-with-sample-counts.xls", cell.type.by.accession.with.sample.counts, sep="\t", row.names=FALSE, col.names=TRUE)


write.table(file = "cell-type-by-accession.xls", cell.type.by.accession, sep="\t", row.names=FALSE, col.names=TRUE)



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
