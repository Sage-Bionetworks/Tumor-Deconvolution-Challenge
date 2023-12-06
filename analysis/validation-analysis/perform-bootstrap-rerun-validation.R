suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(gridExtra))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))
suppressPackageStartupMessages(p_load("reshape2"))

source("../utils.R")

set.seed(1234)

figs.dir <- "figs/"
dir.create(figs.dir, showWarnings = FALSE)

synLogin()

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Read in the rerun predicitons (i.e., where the coarse- and fine-grained datasets are the same)
synId <- "syn22320329"
obj <- synGet(synId, downloadFile=TRUE)
res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

## Read in the bootstraps.rds file
synId <- "syn22344963"
obj <- synGet(synId, downloadFile=TRUE)
bootstraps <- readRDS(obj$path)

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
method.name.col <- "method.name"
round.col <- "submission"

# this defines methods.to.exclude
source("methods-to-exclude.R")

flag <- res.all[,method.name.col] %in% methods.to.exclude
cat(paste0("Excluding methods: ", paste(unique(res.all[flag, method.name.col]), collapse = ", "), "\n"))
res.all <- res.all[!flag,]

## Ensure we have a prediction (even if it is NA) for all cell types in all datasets by all methods
tmp <- unique(res.all[, !(colnames(res.all) %in% c(cell.type.col, prediction.col, measured.col))])
cell.types.by.sub.challenge <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col)])
tmp <- merge(tmp, cell.types.by.sub.challenge, all = TRUE)
measured.tbl <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, measured.col)])
prediction.tbl <- res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, method.name.col, round.col, prediction.col)]
tmp <- merge(tmp, measured.tbl, all = TRUE)
tmp <- merge(tmp, prediction.tbl, all = TRUE)
res.all <- tmp
flag <- !is.na(res.all[, measured.col])
res.all <- res.all[flag, ]

res.all[, cell.type.col] <- as.character(res.all[, cell.type.col])
res.all <- rename.cell.types(res.all, from.col = cell.type.col, to.col = cell.type.col)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))

source("validation-analysis-utils.R")

results <- list()
## for(round in c("1", "2", "3", "latest")) {
for(round in c("1", "2", "3", "latest")) {
# for(round in c("1")) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))

    if(FALSE) {
        res.input <- res.all
        round.col <- "submission"
        round <- "2"
        postfix <- paste0("-round-", round)
    }
    
    results[[round]] <- do.bootstrap.analysis(res.all, bootstraps, 
                                              method.name.col,
                                              subchallenge.col, measured.col, cell.type.col,
                                              dataset.name.col, sample.id.col, prediction.col,
                                              round.col = "submission", round = round)
    
    cat("Bayes factor resultes K <= 3 (or K <= 5) suggests a tie\n")
    for(sub.challenge in sub.challenges) {
        best.team <- results[[round]][["top.performers"]][[sub.challenge]]
        cat(paste0("Top performer for Round ", round, " in ", sub.challenge, " sub-challenge: ",
                   best.team, "\n"))
        tbl <- results[[round]][["bayes.factors"]][[sub.challenge]]
        colnames(tbl)[1] <- "team"
        print(tbl)
        file <- paste0("rerun-validation-bayes-factors-round-", round, "-sc-", sub.challenge, "-vs-",
                       make.names(best.team), ".tsv")
        write.table(file = file, tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    }
}

## save.image(".Rdata.plot.bootstrap")

## Store the above results to Synapse
rounds <- names(results)
    
parent.id <- "syn22320184"

url.base <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
this.script <- "perform-bootstrap-rerun-validation-analysis.R"
script_url <- paste0(url.base, "/", "validation-analysis", "/", this.script)

file <- "bootstrap-validation-results.rds"
saveRDS(results, file)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, executed = script_url, synapseStore = TRUE)
synStore(f)

cat("Exiting successfully\n")
q(status=0)


