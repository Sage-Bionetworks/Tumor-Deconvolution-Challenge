suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(gridExtra))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

source("../utils.R")

set.seed(1234)

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

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
method.name.col <- "method.name"
round.col <- "submission"

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

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")
correlation.methods <- list("pearson" = "pearson", "spearman" = "spearman")

cors <-
    llply(correlation.methods,
          .fun = function(correlation.method) {
              cor.dataset.celltype <-
                  ddply(res.all,
                        .variables = c(subchallenge.col, round.col, cell.type.col, dataset.name.col, method.name.col),
                        .fun = function(df) {
                            val <- NA
                            if(any(is.na(df[, prediction.col]))) {
                                if(!all(is.na(df[, prediction.col]))) {
                                    stop("Some but not all NA\n")
                                }
                            } else {
                                if(var(as.numeric(df[, prediction.col])) == 0) {
                                    val <- 0
                                } else {
                                    val <- cor(df[, measured.col], df[, prediction.col], method = correlation.method)
                                }
                            }
                            data.frame(cor = val)
                        })

              ## Average over dataset
              cor.celltype <-
                  ddply(cor.dataset.celltype,
                        .variables = c(subchallenge.col, round.col, cell.type.col, method.name.col),
                        .fun = function(df) {
                            data.frame(cor = mean(df$cor))
                        })
              cor.celltype
          })

rounds <- c("1", "2", "3", "latest")

## Spot check that the methods have the same scores for coarse- and fine-grained
## NB: some methods may differ between coarse and fine-grained; pick several
## baseline/comparator methods that we know to be the same
check.methods <- c("CIBERSORT", "MCP-counter", "CIBERSORTx")
for(meth in check.methods) {
    coarse.flag <- ( cors[["pearson"]][, method.name.col] == meth ) & ( cors[["pearson"]][, subchallenge.col] == "coarse" )
    meth.res.coarse <- cors[["pearson"]][coarse.flag,]
    fine.flag <- ( cors[["pearson"]][, method.name.col] == meth ) & ( cors[["pearson"]][, subchallenge.col] == "fine" )
    meth.res.fine <- cors[["pearson"]][fine.flag,]
    m <- merge(meth.res.coarse, meth.res.fine, by = c(cell.type.col, round.col))
    m$diff <- m$cor.x - m$cor.y
    m <- m[!is.na(m$cor.x),]
    eps <- 10^-4
    cat(paste0(meth, " max diff between coarse and fine-grained is: ", max(abs(m$diff)), "\n"))
    if(any(abs(m$diff) > eps)) { stop("Max diff exceeded\n") }
}

for(sc in sub.challenges) {
    for(round in rounds) {
        means <- cors[["pearson"]]
        flag <- means[, subchallenge.col] == sc
        means <- means[flag,]

        ## flag <- means[, round.col] == round
        ## means <- means[flag, ]
        
        submitter.tbl <- unique(means[, c(method.name.col, round.col, subchallenge.col), drop = FALSE])
        or <- order(submitter.tbl[, round.col])
        submitter.tbl <- submitter.tbl[or, ]
        submitter.tbl[, round.col] <- as.character(submitter.tbl[, round.col])
        flag <- submitter.tbl[, round.col] %in% c("latest", as.character(round))
        submitter.tbl <- submitter.tbl[flag, ]
        flag <- !duplicated(submitter.tbl[, c(method.name.col, subchallenge.col)], fromLast = FALSE)
        submitter.tbl <- submitter.tbl[flag, ]

        means <- merge(means, submitter.tbl)
        
        g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                                id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        title <- paste0(firstup(sc), "-Grained Sub-Challenge")
        round.text <- ""
        if(round == "latest") {
            round.text <- "Latest Round"
        } else if (round == "1") {
            round.text <- "Round 1"
        } else {
            round.text <- paste0("Latest Round up to Round ", round)
        }
        title <- paste0(title, " (", round.text, ")")
        
        g <- g + ggtitle(title)
        png(paste0("rerun-validation-cell-heatmap-", sc, "-round-", round, ".png"), width = 2 * 480)
        print(g)
        d <- dev.off()
    }
}
