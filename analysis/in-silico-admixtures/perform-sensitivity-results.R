
suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(openxlsx))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(ggpubr))

## sensitivity / spike-in analysis

## Get the sensitivity / spike-in results
synLogin()

synId <- "syn22858611"
validation.results.file <- synGet(synId, downloadFile = TRUE)$path

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
method.name.col <- "method.name"
round.col <- "submission"

res <- fread(validation.results.file, sep=",")
res <- as.data.frame(res)

source("../utils.R")

res <- assign.result.team.names.and.rounds(res, error.fun = warning)

cat("Done assigning team names and rounds\n")

## Read in results from updated CIBERSORTx code that collapses sub-popuations (e.g., mem CD8 T and naive CD8 T)
## into their parental population (e.g., CD8 T) by correctly using sum rather than mean (as was used
## during the challenge). We will apply these changes to the CIBERSORT results (without re-running them).
## These files have columns cell.type, subchallenge, method.name, and revised.orig.ratio.
## method.name = CIBERSORTx (not CIBERSORT!) since CIBERSORTx combines the same sub-populations into parental units
## and was actually re-run. 
## We should revise the old cs result by multiplying it be revised.orig.ratio (i.e., convert mean to sum)
cs.revised.synIds <- list("coarse" = "syn26143304", "fine" = "syn26142832")
## syn26143304: coarse-in-silico-spikeins-csx-all-gene-predictions-revised-orig-ratios.tsv
## syn26142832: fine-in-silico-spikeins-csx-all-gene-predictions-revised-orig-ratios.tsv
## Hmmm ... for some reason, both of these files contain both coarse and fine-grained mappings.
## Ahh ... the different files result from the coarse- vs fine-grained subchallenges / datasets.
## But, CIBERSORTx was run so as to make coarse- or fine-grained predictions for both challenges.
## Proably the files are identical, but to be sure, only take the fine-grained mappings from the
## fine-grained challenge, etc.
revised.dfs <- 
  ldply(cs.revised.synIds, 
        .fun = function(synId) { 
                 obj <- synGet(synId, downloadFile=TRUE)
                 df <- read.table(obj$path, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
               })
colnames(revised.dfs)[1] <- "subchallenge.run"
revised.dfs$method.name <- "CIBERSORT"
revised.dfs <- subset(revised.dfs, subchallenge == subchallenge.run)
revised.dfs <- revised.dfs[, !(colnames(revised.dfs) == "subchallenge.run")]

print(head(res))
print(head(revised.dfs))

orig.nrows <- nrow(res)
res <- merge(as.data.frame(res), as.data.frame(revised.dfs), all.x = TRUE)
new.nrows <- nrow(res)
if(orig.nrows != new.nrows) { stop(paste0("Dimension changed from ", orig.nrows, " rows to ", new.nrows, " rows\n")) }
flag <- !is.na(res$revised.orig.ratio)
if(any(flag)) {
  res[flag, "prediction"] <- res[flag, "prediction"] * res[flag, "revised.orig.ratio"]
}


include.csx.results <- TRUE
if(include.csx.results) {
    cat("Including CIBERSORTx\n")
    ## synIds <- list("coarse" = "syn22725972", "fine" = "syn22726111")
    synIds <- list("coarse" = "syn22340351", "fine" = "syn22340322")
    nms <- names(synIds)
    names(nms) <- nms
    csx.res <-
        ldply(nms,
              .fun = function(nm) {
                  synId <- synIds[[nm]]
                  csx.results.file <- synGet(synId, downloadFile = TRUE)$path
                  tbl <- read.table(csx.results.file, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
                  flag <- tbl$subchallenge == nm
                  tbl <- tbl[flag, ]
                  tbl$objectId <- "CIBERSORTx"
                  tbl$repo_name <- "CIBERSORTx"
                  tbl$comparator <- TRUE
                  tbl$submitterId <- NA
                  tbl1 <- tbl
                  tbl2 <- tbl
                  tbl1$submission <- "1"
                  tbl2$submission <- "latest"                  
                  rbind(tbl1, tbl2, stringsAsFactors = FALSE)
              })
    measured <- unique(res[, c("subchallenge", "dataset.name", "sample.id", "cell.type", "measured")])
    csx.res <- merge(csx.res, measured)
  
    cols <- intersect(colnames(res), colnames(csx.res))

    res <- rbind(res[, cols], csx.res[, cols])
}


flag <- res$sample.id == "Breast"
res[flag, "sample.id"] <- "BRCA"

cat("Reading in spike-in annotations\n")

## Read in the spike-in annotations (i.e., indicating what the spike-in population was)
synIds <- list("coarse" = "syn21763908", "fine" = "syn21763907")
spikein.annotation <-
    ldply(synIds,
          .fun = function(synId) {
              admixture.file <- synGet(synId, downloadFile = TRUE)$path
              tbl <- read.table(admixture.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
              tbl <- subset(tbl, spike.in.pop == sample)
              ret <- unique(tbl[, c("dataset.name", "sample.id", "spike.in.pop")])
              colnames(ret) <- c(dataset.name.col, sample.id.col, cell.type.col)
              return(ret)
              ## For coarse-grained, need to collapse the specific cell types comprising the parental population
              ## e.g., myeloid lineage = DCs + macrophages + myeloid dendritic cells
              ddply(tbl,
                    .variables = c("dataset.name", "sample.id", "spike.in.pop"),
                    .fun = function(df) {
                        ## Need this rounding because sometimes the sum is like 0.00499999999999999
                        sm <- round(sum(df$measured), digits = 5)
                        data.frame(measured = sm)
                    })
          })
colnames(spikein.annotation)[1] <- subchallenge.col

cat("Merging spike-in annotations\n")

## We only care about the results for the spike in populations, so merge to effectively limit to those
res <- merge(res, spikein.annotation)


## Let's exclude memory.B.cells (which always have measured == 0, which causes problems with correlation)
res <- subset(res, !(cell.type == "memory.B.cells"))

cat("Renaming cell types\n")

res[, cell.type.col] <- as.character(res[, cell.type.col])
res <- rename.cell.types(res, from.col = cell.type.col, to.col = cell.type.col)

## We can have redudant cell types if the underlying samples are different, e.g.,
## monocytic.lineage = Macrophages
## monocytic.lineage = Monocytes
## monocytic.lineage = DCs
## we want their measured sum (which corresponds to the 3 samples in this example).
## The predictions should all be the same (i.e., for monocytic.lineage)
print(colnames(res))

cat("Summing coarse-grained cell types\n")

all.res.cols <- c(subchallenge.col, dataset.name.col, sample.id.col, cell.type.col, method.name.col, round.col, "comparator")
res <- res[, c(all.res.cols, "measured", "prediction")]
res <-
    ddply(res,
          .variables = all.res.cols,
          .fun = function(df) {
              if((nrow(df) > 1) && (df[1, subchallenge.col] == "fine")) {
                  print(df)
                  stop("Was not expecting multiple repeats for fine grained challenge\n")
              }
              preds <- df[, "prediction"]
              if(!all(is.na(preds)) && !all(preds == preds[1])) {
                  print(df)
                  stop("Was not expecting different predicted results\n")
              }
              ## Need this rounding because sometimes the sum is like 0.00499999999999999              
              measured.sum <- round(sum(df$measured), digits = 5)
              if(measured.sum > 1) {
                  print(df)
                  stop("df > 1")
              }
              data.frame(measured = measured.sum, prediction = preds[1])
          })

cat("Merging in metadata\n")

## Add the mixture.type annotation
annotations <- get.in.silico.metadata()
annotations <- annotations[, c("dataset.name", "mixture.type", "subchallenge")]

res <- safe.merge(res, annotations)

sig.pval.cutoff <- 0.01
plot.pval.cutoff <- sig.pval.cutoff

perform.sensitivity.analysis <- function(res.input, postfix,
                                         method.name.col, 
                                         subchallenge.col, 
                                         round.col, round = "latest",
                                         mixture.col = "mixture.type",
                                         cell.type.col = "cell.type",
                                         dataset.name.col = "dataset.name") {

    submitter.tbl <- unique(res.input[, c(method.name.col, round.col, subchallenge.col), drop = FALSE])
    or <- order(submitter.tbl[, round.col])
    submitter.tbl <- submitter.tbl[or, ]
    submitter.tbl[, round.col] <- as.character(submitter.tbl[, round.col])
    flag <- submitter.tbl[, round.col] %in% c("latest", as.character(round))
    submitter.tbl <- submitter.tbl[flag, ]
    flag <- !duplicated(submitter.tbl[, c(method.name.col, subchallenge.col)], fromLast = FALSE)
    submitter.tbl <- submitter.tbl[flag, ]

    res <- merge(res.input, submitter.tbl, by = c(method.name.col, subchallenge.col, round.col))

    round.text <- ""
    if(round == "latest") {
        round.text <- "Latest Round"
    } else if (round == "1") {
        round.text <- "Round 1"
    } else {
        round.text <- paste0("Latest Round up to Round ", round)
    }

    title.postfix <- round.text

    res.matrices <-
        dlply(res, .variables = c(subchallenge.col),
              .fun = function(res.sc) {
                  ddply(res.sc,
                        .variables = c(cell.type.col, dataset.name.col, method.name.col, mixture.col),
                        .parallel = FALSE,
                        .fun = function(df) {
                            if(all(is.na(df$prediction))) {
                                return(data.frame(min.diff.prop = NA))
                            }
                            if(any(is.na(df$prediction))) {
                                print(df)
                                stop("Was not expecting some NA and some not\n")
                            }

                            df <- df[order(df$measured),,drop=F]
                            sci <- formatC(as.numeric(as.character(df$measured)), format="e", digits=2)
                            df$measured <- factor(sci, levels = unique(sci))
                            
                            cmps <- as.data.frame(compare_means(prediction ~ measured,  data = df))
                            cmps <- cmps[cmps$group1 == "0.00e+00",]
                            ret <- cmps[, c("group1", "group2", "p")]
                            min.diff.prop <- 1
                            if(any(ret$p < sig.pval.cutoff)) {
                                use.first.sig.diff <- FALSE
                                if(use.first.sig.diff) {
                                    ## Take the lowest spike-in level that is significantly different
                                    flag <- ret$p < sig.pval.cutoff
                                    min.diff.prop <- min(as.numeric(ret[flag,"group2"]))
                                } else {
                                    ## Take the lowest spike-in level such that it and all larger levels are
                                    ## significantly different
                                    ret <- ret[order(as.numeric(ret$group2), decreasing=TRUE),]
                                    for(idx in 1:nrow(ret)) {
                                        if(all(ret[1:idx, "p"] < sig.pval.cutoff)) {
                                            min.diff.prop <- as.numeric(ret[idx, "group2"])
                                        } else {
                                            break
                                        }
                                    }
                                }
                                if(min.diff.prop > 1) {
                                    print(df)
                                    print(ret)
                                    stop("min.diff.prop > 1")
                                }
                            }
                            data.frame(min.diff.prop = min.diff.prop)
                        })
              })

    res.list <-
        list("res" = res,
             "res.matrices" = res.matrices)
    
    return(res.list)
}

results <- list()
rounds <- c("1", "2", "3", "latest")
for(round in rounds) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))

    results[[round]] <- perform.sensitivity.analysis(res, postfix = postfix, method.name.col = method.name.col,
                                                     subchallenge.col = subchallenge.col,
                                                     round.col = "submission", round = round,
                                                     mixture.col = "mixture.type", cell.type.col = cell.type.col)
                                         
}

cat("Saving image\n")
save.image(".Rdata.sensitivity")
cat("Done saving image\n")

parent.id <- "syn22951635"

url.base <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
this.script <- "perform-sensitivity-analysis.R"
script_url <- paste0(url.base, "/", "in-silico-admixtures", "/", this.script)

file <- "sensitivity-analysis-results.rds"
saveRDS(results, file)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, executed = script_url, synapseStore = TRUE)
synStore(f)

file <- "sensitivity-analysis-data.rds"
saveRDS(res, file)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, executed = script_url, synapseStore = TRUE)
synStore(f)

cat("Exiting successfully\n")
q(status = 0)


## working

my.format <- function(x) {
    ifelse(x == 0, 0,
           ifelse(x < 0.01, formatC(x, format="e", digits=1),
                  formatC(x, format="f", digits=2, drop0trailing = TRUE)))
}
    
res.matrices <- results[[1]]$res.matrices
nms <- names(res.matrices)
names(nms) <- nms
all.plots <- 
    llply(nms,
          .fun = function(nm) {
              res <- res.matrices[[nm]]
              subplots <-
                  dlply(res, .variables = c("mixture.type"),
                        .fun = function(res.ds) {

                            print(head(res.ds))
                            df <- acast(res.ds, as.formula(paste0(method.name.col, " ~ ", cell.type.col)), value.var = "min.diff.prop", fill = NA)
                            res.ds <- reshape2::melt(as.matrix(df))
                            colnames(res.ds) <- c(method.name.col, cell.type.col, "min.diff.prop")
                            print(head(res.ds))

                            id.var <- method.name.col
                            
                            res.ds$label <- as.numeric(as.character(100 * res.ds$min.diff.prop))


                            
                            ## Cast from long to matrix form to introduce NAs for missing values,
                            ## then melt back to long form
                            
                            g <- plot.cell.type.correlation.heatmap(res.ds, show.corr.text = TRUE,
                                                                    id.var = method.name.col, cell.type.var = cell.type.col,
                                                                    cor.var = "label", formatter = my.format,
                                                                    col.summary.fun = "min", cor.type.label = "LoD (Percent)",
                                                                    limits = c(0, 100),
                                                                    order.decreasing = TRUE)
                            ## limits = c(log2(min(res$measured[res$measured > 0])), log2(1))
                            print(min(res.ds$label))
                            print(max(res.ds$label))
                            limits = c(log2(10^-4), log2(1))
                            my.labeller <- function(x) { my.format(100*(2^x)) }
                            g <- g + scale_fill_gradient2("LoD (Percent)", labels = my.labeller, limits = limits, 
                                                          low = "red", high = "blue", mid = "white", na.value = "black")
                            
                            g
                        })
          })


                            res.ds$cell.type <- factor(res.ds$cell.type, levels = cell.type.levels)
                            res.ds[, id.var] <- factor(res.ds[, id.var], levels = id.levels)
                            
                            g <- ggplot(data = res.ds, aes(y = repo_name, x = cell.type, fill = log2(min.diff.prop)))
                            g <- g + geom_tile()
                            g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                           text = element_text(size = 20),
                                           title = element_text(size = 20))                                           
                            g <- g + ylab("Method") + xlab("")
                            limits = c(min(unname(log2(spike.ins[spike.ins>0]))),log2(1))
                            my.labeller <- function(x) { my.format(100*(2^x)) }



nms <- names(res.matrices)
names(nms) <- nms
all.plots <- 
    llply(nms,
          .fun = function(nm) {
              res <- res.matrices[[nm]]
              subplots <-
                  dlply(res, .variables = c("mixture.type"),
                        .fun = function(res.ds) {

                            var <- "min.diff.prop"
                            id.var <- "repo_name"

                            res.ds$cell.type <- as.character(res.ds$cell.type)
                            res.ds[,id.var] <- as.character(res.ds[,id.var])

                            ## Add mean row and column
                            row.summary.name <- "mean"
                            row.summary.fun <- mean
                            col.summary.name <- "min"
                            col.summary.fun <- min
                            cell.type.summaries <- ddply(res.ds, .variables = c("cell.type"),
                                                         .fun = function(tmp) {
                                                             ret <- data.frame(id = col.summary.name, cor = col.summary.fun(tmp[, var], na.rm=TRUE))
                                                             colnames(ret)[1] <- id.var
                                                             colnames(ret)[2] <- var                                 
                                                             ret
                                                         })
                            cell.type.summaries <- cell.type.summaries[order(cell.type.summaries[, var], decreasing = TRUE),]
    
                            method.summaries <- ddply(res.ds, .variables = c(id.var),
                                                      .fun = function(tmp) {
                                                          ret <- data.frame(cell.type = row.summary.name, cor = row.summary.fun(tmp[, var], na.rm=TRUE))
                                                          colnames(ret)[2] <- var
                                                          ret
                                                      })
    
                            method.summaries <- method.summaries[order(method.summaries[, var], decreasing = TRUE),]


                            cell.type.levels <- c(cell.type.summaries$cell.type[cell.type.summaries$cell.type != row.summary.name], row.summary.name)
                            id.levels <- c(col.summary.name, method.summaries[method.summaries[, id.var] != col.summary.name, id.var])

                            
                            res.ds <- rbind(res.ds[, c(id.var, "cell.type", var)], cell.type.summaries[, c(id.var, "cell.type", var)])
                            res.ds <- rbind(res.ds, method.summaries[, c(id.var, "cell.type", var)])    

                            
                            ## res.ds$label <- formatC(as.numeric(as.character(res.ds$min.diff.prop)), format="e", digits=1)
                            ## res.ds$label <- formatC(as.numeric(as.character(100 * res.ds$min.diff.prop)), format="e", digits=1)

                            res.ds$label <- my.format(as.numeric(as.character(100 * res.ds$min.diff.prop)))
                            res.ds$cell.type <- factor(res.ds$cell.type, levels = cell.type.levels)
                            res.ds[, id.var] <- factor(res.ds[, id.var], levels = id.levels)
                            
                            g <- ggplot(data = res.ds, aes(y = repo_name, x = cell.type, fill = log2(min.diff.prop)))
                            g <- g + geom_tile()
                            g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                           text = element_text(size = 20),
                                           title = element_text(size = 20))                                           
                            g <- g + ylab("Method") + xlab("")
                            limits = c(min(unname(log2(spike.ins[spike.ins>0]))),log2(1))
                            my.labeller <- function(x) { my.format(100*(2^x)) }
                            g <- g + scale_fill_gradient2("LoD (Percent)", labels = my.labeller, limits = limits, 
                                                          low = "red", high = "blue", mid = "white", na.value = "black")
                            g <- g + geom_text(aes(label = label))
                            mt <- as.character(res.ds$mixture.type[1])
                            title <- paste0(mt, " Admixtures")
                            g <- g + ggtitle(title)
                        })
              subplot.names <- names(subplots)
              names(subplot.names) <- subplot.names
              l_ply(subplot.names,
                    .fun = function(sb.nm) {
                        mt <- sb.nm
                        prefix <- mt
                        title <- paste0(mt, " Admixtures")
                        g <- subplots[[sb.nm]]
                        g <- g + ggtitle(paste0(firstup(nm), "-Grained Sub-Challenge (", title, ")"))
                        g <- g + theme(text = element_text(size = 20), title = element_text(size = 20))
                        file <- paste0("validation-lod-summary-", prefix, "-", nm, ".png")
                        png(file, width = 2 * 480)
                        print(g)
                        d <- dev.off()
                    })
              ## file <- paste0("validation-lod-summary-", nm, ".png")
              ## png(file, width = 2 * 480)
              ## title <- paste0(firstup(nm), "-Grained Sub-Challenge (Validation)")
              ## g <- do.call("grid.arrange", c(subplots, nrow = 2, top = textGrob(title, gp=gpar(fontsize=25))))
              ## grid.draw(g)
              ## d <- dev.off()
              subplots
          })

g1 <- (all.plots[["coarse"]][["Random"]])
g2 <- (all.plots[["fine"]][["Random"]])
g1 <- g1 + ggtitle("Coarse-Grained Sub-Challenge")
g1 <- g1 + theme(text = element_text(size = 20), title = element_text(size = 20))
g2 <- g2 + ggtitle("Fine-Grained Sub-Challenge")
g2 <- g2 + theme(text = element_text(size = 20), title = element_text(size = 20))
file <- paste0("validation-lod-summary-random", ".png")
png(file, width = 2 * 480, height = 2 * 480)
title <- "Random Admixtures"
## g <- do.call("grid.arrange", c(list(g1,g2), nrow = 2, top = title))
g <- grid.arrange(g1, g2, nrow = 2, top = textGrob(title, gp=gpar(fontsize=25)))
grid.draw(g)
d <- dev.off()

g1 <- (all.plots[["coarse"]][["Random"]])
g2 <- (all.plots[["coarse"]][["Biological"]])
g1 <- g1 + ggtitle("Random Admixtures")
g1 <- g1 + theme(text = element_text(size = 20), title = element_text(size = 20))
g2 <- g2 + ggtitle("Biological Admixtures")
g2 <- g2 + theme(text = element_text(size = 20), title = element_text(size = 20))
file <- paste0("validation-lod-summary-coarse", ".png")
png(file, width = 2 * 480, height = 2 * 480)
title <- "Coarse-Grained Sub-Challenge"
## g <- do.call("grid.arrange", c(list(g1,g2), nrow = 2, top = title))
g <- grid.arrange(g1, g2, nrow = 2, top = textGrob(title, gp=gpar(fontsize=25)))
grid.draw(g)
d <- dev.off()

g1 <- (all.plots[["fine"]][["Random"]])
g2 <- (all.plots[["fine"]][["Biological"]])
g1 <- g1 + ggtitle("Random Admixtures")
g1 <- g1 + theme(text = element_text(size = 20), title = element_text(size = 20))
g2 <- g2 + ggtitle("Biological Admixtures")
g2 <- g2 + theme(text = element_text(size = 20), title = element_text(size = 20))
file <- paste0("validation-lod-summary-fine", ".png")
png(file, width = 2 * 480, height = 2 * 480)
title <- "Fine-Grained Sub-Challenge"
## g <- do.call("grid.arrange", c(list(g1,g2), nrow = 2, top = title))
g <- grid.arrange(g1, g2, nrow = 2, top = textGrob(title, gp=gpar(fontsize=25)))
grid.draw(g)
d <- dev.off()


nms <- names(spike.in.res)
names(nms) <- nms
l_ply(nms,
      .fun = function(nm) {
          res <- spike.in.res[[nm]]
          d_ply(res,
                .variables = c("cell.type", "repo_name"),
                .fun = function(df) {
                    ct <- as.character(df$cell.type[1])
                    subplots <-
                        dlply(df,
                              .variables = c("mixture.type"),
                              .fun = function(df.ds) {


                                  df.ds <- df.ds[order(df.ds$measured),]
                                  df.ds$val <- as.numeric(as.character(df.ds$measured)) * 100
                                  ## sci <- formatC(as.numeric(as.character(df.ds$measured))*100, format="e", digits=2)
                                  sci <- my.format(df.ds$val)
                                  if(!use.segmented) {
                                      df.ds$val <- factor(sci, levels = unique(sci))
                                      cmps <- as.data.frame(compare_means(prediction ~ val,  data = df.ds))
                                      ## cmps <- cmps[cmps$group1 == "0.00e+00",]
                                      cmps <- cmps[cmps$group1 == "0",]
                                      display.p.cutoff <- 0.1
                                      display.p.cutoff <- 0.01
                                      cmps <- cmps[cmps$p < display.p.cutoff,]
                                      my_comparisons <-
                                          llply(1:nrow(cmps),
                                                .fun = function(i) c(cmps[i,2], cmps[i,3]))
                                  }
                                  
                                  g <- ggplot(data = df.ds, aes(x = val, y = prediction))
                                  g <- g + geom_boxplot()
                                  if(!use.segmented) {
                                      g <- g + stat_compare_means(comparisons = my_comparisons, size = 5)
                                  }
                                  g <- g + geom_beeswarm()
                                  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                                  xlab <- paste0(fix.string(ct), " Spike-in (%)")
                                  g <- g + xlab(xlab) + ylab("Score")
                                  mt <- as.character(df.ds$mixture.type[1])
                                  title <- paste0(mt, " Admixtures")
                                  g <- g + ggtitle(title)
                                  g
                              })
                    meth <- as.character(df$repo_name[1])
                    subplot.names <- names(subplots)
                    names(subplot.names) <- subplot.names
                    l_ply(subplot.names,
                          .fun = function(sb.nm) {
                              mt <- sb.nm
                              prefix <- mt
                              title <- paste0(mt, " Admixtures")
                              g <- subplots[[sb.nm]]
                              g <- g + ggtitle(paste0(meth, ": ", ct, "\n", firstup(nm), "-Grained Sub-Challenge (", title, ")"))
                              g <- g + theme(text = element_text(size = 20), title = element_text(size = 20))
                              file <- paste0("validation-lod-ct-", prefix, "-", ct, "-", meth, "-", nm)
                              file <- gsub(file, pattern = "\\.", replacement = "-")
                              file <- paste0(file, ".png")
                              png(file, width = 2 * 480)
                              print(g)
                              d <- dev.off()
                    })

                    top <- paste0(meth, ": ", ct, "\n", firstup(nm), "-Grained Sub-Challenge")
                    plt <- do.call("grid.arrange", c(subplots, nrow = 2, top = top))
                    file <- make.names(paste0("validation-lod-ct-", ct, "-", meth, "-", nm))
                    file <- gsub(file, pattern = "\\.", replacement = "-")
                    file <- paste0(file, ".png")
                    print(file)
                    png(file, width = 2 * 480)
                    grid.draw(plt)
                    d <- dev.off()
                })
      })

cat("Exiting successfully\n")
q(status=0)

plot.cell.type.score.heatmap <- function(df, score.col = "prediction",
                                         normalized.score = FALSE, nrow = 1, method.name.col = "method.name") {

    ## This re-ordering o works with as.table = FALSE to respect the
    ## the ordering we want from cell.type.levels
    if(nrow > 1) {
        method.levels <- sort(unique(df[, method.name.col]), decreasing = TRUE)
        ntot <- length(method.levels)
        ncol <- ceiling(ntot / nrow)
        o <- unlist(llply(1:(nrow-1), .fun = function(i) (i*ncol):(1+((i-1)*ncol))))
        o <- c(o, ntot:(max(o)+1))
        df[, method.name.col] <- factor(df[, method.name.col], levels = method.levels[o])
    }
    g <- ggplot(data = df, aes_string(y = "sample.id", x = "cell.type", fill = score.col))
    g <- g + geom_tile()
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   text = element_text(size = 16))
##                   title = element_text(size = 8))
    g <- g + ylab("Purified Sample") + xlab("Predicted Cell Type")
    if(normalized.score) {
        g <- g + scale_fill_gradient2("Normalized\nPrediction", 
                                      low = "red", high = "blue", mid = "white", na.value = "black")
    } else {
        g <- g + scale_fill_gradient2("Prediction", limits = c(0,1),
                                      low = "red", high = "blue", mid = "white", na.value = "black")
    }
    g <- g + facet_wrap(method.name.col, as.table = FALSE, nrow = nrow)
    g
}

