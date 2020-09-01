
suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggpubr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(ggbeeswarm))

sig.pval.cutoff <- 0.01
plot.pval.cutoff <- sig.pval.cutoff

synLogin()
synId <- "syn21739521"
cat(paste0("NB: using Version 2 here! Don't know how these versions got screwed up!\n"))
in.silico.results.file <- synGet(synId, downloadFile = TRUE, version = 2)$path
res <- read.table(in.silico.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

## Read in the cibersortX results
## synId <- "syn22331779"
## csx.in.silico.results.file <- synGet(synId, downloadFile = TRUE)$path
## csx.res <- read.table(csx.in.silico.results.file, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

## Read in the fine-grained cibersortX results
synId <- "syn22340322"
fine.csx.in.silico.results.file <- synGet(synId, downloadFile = TRUE)$path
fine.csx.res <- read.table(fine.csx.in.silico.results.file, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

## Read in the coarse-grained cibersortX results
synId <- "syn22340351"
coarse.csx.in.silico.results.file <- synGet(synId, downloadFile = TRUE)$path
coarse.csx.res <- read.table(coarse.csx.in.silico.results.file, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

csx.res <- rbind(fine.csx.res, coarse.csx.res)

source("../utils.R")

val.metadata <- get.in.silico.metadata()
res <- merge(res, val.metadata, all.x = TRUE, by = c("dataset.name", "subchallenge"))

## We want to use the datasets that did not have cancer
res <- res[is.na(res$tumor.type),]

csx.res <- subset(csx.res, dataset.name %in% res$dataset.name)

merge.measured.manually <- TRUE
if(merge.measured.manually) {
    synId <- "syn21763908"
    coarse.gt.file <- synGet(synId, downloadFile = TRUE)$path
    coarse.gt <- read.table(coarse.gt.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
    coarse.gt$subchallenge <- "coarse"

    synId <- "syn21763907"
    fine.gt.file <- synGet(synId, downloadFile = TRUE)$path
    fine.gt <- read.table(fine.gt.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
    fine.gt$subchallenge <- "fine"

    gt <- rbind(coarse.gt, fine.gt)
    gt <- gt[, c("dataset.name", "sample.id", "subchallenge", "sample", "measured")]
    colnames(gt) <- c("dataset.name", "sample.id", "subchallenge", "cell.type", "measured")
    csx.res <- merge(csx.res, unique(gt))


} else {
    ## Merge the measured/ground truth into CIBERSORTx predictions from those
    ## measured values already included in the predictions for the
    ## non-CIBERSORTx methods
    csx.res <- merge(csx.res, unique(res[, c("dataset.name", "sample.id", "subchallenge", "mixture.type", "cell.type", "measured")]))
}

csx.res <- merge(csx.res, val.metadata, all.x = TRUE, by = c("dataset.name", "subchallenge"))
csx.res$repo_name <- "CIBERSORTx"

res <- csx.res

if(FALSE) {
cols <- intersect(colnames(csx.res), colnames(res))
res <- rbind(res[, cols], csx.res[, cols])
}

## This would effectively limit to the same datasets:
## res <- subset(res, dataset.name %in% c("A", "B", "G", "H"))

## synId <- "syn21752552"
synId <- "syn21763908"
in.silico.coarse.admixture.file <- synGet(synId, downloadFile = TRUE)$path
in.silico.coarse.admixtures <- read.table(in.silico.coarse.admixture.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

## synId <- "syn21752551"
synId <- "syn21763907"
in.silico.fine.admixture.file <- synGet(synId, downloadFile = TRUE)$path
in.silico.fine.admixtures <- read.table(in.silico.fine.admixture.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

## Plot scores for cell type x as a function of spike in for cell type x for each cell type, method, subchallenge

trans <-
    list("baseline_method1" = "CIBERSORT",
         "baseline_method2" = "MCP-counter",
         "baseline_method3" = "quanTIseq",
         "baseline_method4" = "xCell",
         "baseline_method5" = "EPIC",
         "baseline_method6" = "TIMER",
         "CelEst" = "EPIC (applicant)")
for(nm in names(trans)) {
    flag <- grepl(res$repo_name, pattern = nm)
    res[flag, "repo_name"] <- trans[[nm]]
}

spike.ins <- c(1/(2^seq(from=2,to=12,by=1)), 0)
names(spike.ins) <- spike.ins

## Begin define map
synLogin()

## If we don't have a reliable annotation of which cell type is the
## spike-in within a given admixture, we have to guess it from the fact
## that it has one of the spike-in values. This sounds simple, but is
## complicated by the fact that, for coarse grained, we spiked in a coarse
## popuation, which has multiple constituent populations
infer.spike.in.cell.types <- FALSE

if(infer.spike.in.cell.types) {
    ## Load the TPM validation data
    ##synId <- "syn21574299"
    ##obj <- synGet(synId, downloadFile = TRUE)
    ##cpm.expr <- read.table(obj$path, sep = "\t", header = TRUE)
    
    synId <- "syn21576632"
    obj <- synGet(synId, downloadFile = TRUE)
    cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)
    
    samples <- colnames(cpm.expr)
    ## Exclude the admixtures samples
    flag <- grepl(samples, pattern="RM")
    samples <- samples[!flag]
    flag <- grepl(samples, pattern="BM")
    samples <- samples[!flag]
    samples <- samples[!(samples == "Gene")]
    names(samples) <- samples
    
    ## Map populations to samples
    populations <- samples
    populations <- gsub(populations, pattern="_1", replacement="")
    populations <- gsub(populations, pattern="_2", replacement="")
    pop.df <- data.frame(population = unname(populations), sample = names(populations))
    
    ## Translate the population names to those we use in the challenge
    fine.grained.map <- list(
        "Breast" = "Breast",
        "CRC" = "CRC",
        "Dendritic_cells" = "myeloid.dendritic.cells",
        "Endothelial_cells" = "endothelial.cells",
        "Fibroblasts" = "fibroblasts",
        "Macrophages" = "macrophages",
        "Memory_CD4_T_cells" = "memory.CD4.T.cells",
        "Memory_CD8_T_cells" = "memory.CD8.T.cells",
        "Monocytes" = "monocytes",
        "NK_cells" = "NK.cells",
        "Naive_B_cells" = "naive.B.cells",
        "Naive_CD4_T_cells" = "naive.CD4.T.cells",
        "Naive_CD8_T_cells" = "naive.CD8.T.cells",
        "Neutrophils" = "neutrophils",
        "Tregs" = "regulatory.T.cells")
    
    coarse.grained.map <- list(
        "Breast" = "Breast",
        "CRC" = "CRC",
        "Dendritic_cells" = "monocytic.lineage",
        "Endothelial_cells" = "endothelial.cells",
        "Fibroblasts" = "fibroblasts",
        "Macrophages" = "monocytic.lineage",
        "Memory_CD4_T_cells" = "CD4.T.cells",
        "Memory_CD8_T_cells" = "CD8.T.cells",
        "Monocytes" = "monocytic.lineage",
        "NK_cells" = "NK.cells",
        "Naive_B_cells" = "B.cells",
        "Naive_CD4_T_cells" = "CD4.T.cells",
        "Naive_CD8_T_cells" = "CD8.T.cells",
        "Neutrophils" = "neutrophils",
        "Tregs" = "CD4.T.cells")
    
    
    fine.grained.map.df <- data.frame(population = names(fine.grained.map), challenge.population = unname(unlist(fine.grained.map)))
    fine.grained.pop.df <- merge(pop.df, fine.grained.map.df)
    
    coarse.grained.map.df <- data.frame(population = names(coarse.grained.map), challenge.population = unname(unlist(coarse.grained.map)))
    coarse.grained.pop.df <- merge(pop.df, coarse.grained.map.df)
    
    ## Separate the samples into two batches
    fine.grained.pop.df <- fine.grained.pop.df[order(fine.grained.pop.df$sample),]
    fine.grained.pop1.df <- ddply(fine.grained.pop.df, .variables = c("population"), .fun = function(df) df[1,,drop=F])
    fine.grained.pop2.df <- ddply(fine.grained.pop.df, .variables = c("population"), .fun = function(df) df[min(2,nrow(df)),,drop=F])
    
    coarse.grained.pop.df <- coarse.grained.pop.df[order(coarse.grained.pop.df$sample),]
    coarse.grained.pop1.df <- ddply(coarse.grained.pop.df, .variables = c("population"), .fun = function(df) df[1,,drop=F])
    coarse.grained.pop2.df <- ddply(coarse.grained.pop.df, .variables = c("population"), .fun = function(df) df[min(2,nrow(df)),,drop=F])
} ## if(infer.spike.in.cell.types)
## End define map

spike.in.res <- list()
spike.in.res[["coarse"]] <- subset(res, subchallenge == "coarse")
## We can have redudant cell types if the underlying samples are different, e.g.,
## monocytic.lineage = Macrophages
## monocytic.lineage = Monocytes
## monocytic.lineage = DCs
## we want their measured sum (which corresponds to the 3 samples in this example).
## The predictions should all be the same (i.e., for monocytic.lineage)
spike.in.res[["coarse"]] <-
    ddply(spike.in.res[["coarse"]],
          .variables = c("dataset.name", "sample.id", "repo_name", "cell.type"),
          .fun = function(df) {
              if(!all(df$prediction == df$prediction[1])) {
                  print(df)
                  stop("Got different measured values\n")
              }
              df$measured <- round(sum(df$measured), digits=5)
              df[1,,drop=F]
          })

if(infer.spike.in.cell.types) {
    spike.in.res[["coarse"]] <- subset(spike.in.res[["coarse"]], measured %in% unname(spike.ins))
    spike.in.res[["coarse"]] <- subset(spike.in.res[["coarse"]], cell.type %in% coarse.grained.map.df$challenge.population)
    if(FALSE) {
        spike.in.res[["coarse"]] <-
            ddply(spike.in.res[["coarse"]],
                  .variables = c("dataset.name", "sample.id"),
                  .fun = function(df) {
                      tmp <- unique(df[, c("cell.type", "measured")])
                      ## Find the cell type that comes closest to a spike.in value
                      si <- expand.grid(cell.type = tmp$cell.type, spike.in = unname(spike.ins))
                      tmp <- merge(tmp, si)
                      tmp$diff <- abs(tmp$measured - tmp$spike.in)
                      spike.in.pop <- as.character(tmp[which.min(tmp$diff), "cell.type"])
                      spike.in.prop <- as.numeric(tmp[which.min(tmp$diff), "spike.in"])
                      df <- df[df$cell.type %in% spike.in.pop,]
                      df$measured <- spike.in.prop
                      df
                  })
    }
} else {
    tmp <- na.omit(unique(in.silico.coarse.admixtures[, c("dataset.name", "sample.id", "spike.in.pop", "spike.in.prop")]))
    spike.in.res[["coarse"]] <- merge(spike.in.res[["coarse"]], tmp, by.x = c("dataset.name", "sample.id", "cell.type"),
                                  by.y = c("dataset.name", "sample.id", "spike.in.pop"))
} ## infer.spike.in.cell.types

spike.in.res[["fine"]] <- subset(res, subchallenge == "fine")
if(infer.spike.in.cell.types) {
    spike.in.res[["fine"]] <- subset(spike.in.res[["fine"]], measured %in% unname(spike.ins))
    spike.in.res[["fine"]] <- subset(spike.in.res[["fine"]], cell.type %in% fine.grained.map.df$challenge.population)
} else {
    tmp <- na.omit(unique(in.silico.fine.admixtures[, c("dataset.name", "sample.id", "spike.in.pop", "spike.in.prop")]))
    spike.in.res[["fine"]] <- merge(spike.in.res[["fine"]], tmp, by.x = c("dataset.name", "sample.id", "cell.type"),
                                  by.y = c("dataset.name", "sample.id", "spike.in.pop"))
} # if(infer.spike.in.cell.types)

res.matrices <-
    llply(spike.in.res,
          .fun = function(res) {
              ddply(res,
                    .variables = c("subchallenge", "cell.type", "dataset.name", "repo_name", "mixture.type"),
                    .fun = function(df) {
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
                        }
                        data.frame(min.diff.prop = min.diff.prop)
                    })
          })

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

nms <- names(res.matrices)
names(nms) <- nms
all.plots <- 
    llply(nms,
          .fun = function(nm) {
              res <- res.matrices[[nm]]
              subplots <-
                  dlply(res, .variables = c("mixture.type"),
                        .fun = function(res.ds) {

                            ## res.ds$label <- formatC(as.numeric(as.character(res.ds$min.diff.prop)), format="e", digits=1)
                            my.format <- function(x) ifelse(x < sig.pval.cutoff, formatC(x, format="e", digits=1),
                                                            formatC(x, format="f", digits=2, drop0trailing = TRUE))
                            ## res.ds$label <- formatC(as.numeric(as.character(100 * res.ds$min.diff.prop)), format="e", digits=1)
                            res.ds$label <- my.format(as.numeric(as.character(100 * res.ds$min.diff.prop)))
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
                        g <- g + ggtitle(paste0(firstup(nm), "-Grained Sub-Challenge (Validation; ", title, ")"))
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
title <- "Random Admixtures (Validation)"
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
title <- "Coarse-Grained Sub-Challenge (Validation)"
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
title <- "Fine-Grained Sub-Challenge (Validation)"
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
                    subplots <-
                        dlply(df,
                              .variables = c("mixture.type"),
                              .fun = function(df.ds) {

                                  df.ds <- df.ds[order(as.numeric(df.ds$measured)),]
                                  ## df.ds$measured <- round(df.ds$measured, digits=3)
                                  use.scientific.notation <- FALSE
                                  if(use.scientific.notation) {
                                      sci <- formatC(as.numeric(as.character(df.ds$measured)), format="e", digits=2)
                                      df.ds$measured <- factor(sci, levels = unique(sci))
                                  } else {
                                      df.ds$measured <- 100 * as.numeric(as.character(df.ds$measured))
                                      df.ds$measured <- factor(df.ds$measured, levels = unique(df.ds$measured))
                                  }
                                  cmps <- as.data.frame(compare_means(prediction ~ measured,  data = df.ds))

                                  if(use.scientific.notation) {                                  
                                      cmps <- cmps[as.character(cmps$group1) == "0.00e+00",]
                                  } else {
                                      cmps <- cmps[as.character(cmps$group1) == "0",]
                                  }
                                  cmps <- cmps[order(as.numeric(as.character(cmps$group2)), decreasing=FALSE),]
                                  cmps$indx <- 1:nrow(cmps)
                                  sig.cmps <- cmps[cmps$p < plot.pval.cutoff,]
                                  my_comparisons <-
                                      llply(1:nrow(sig.cmps),
                                            .fun = function(i) c(sig.cmps[i,2], sig.cmps[i,3]))

                                  ## There are a lot of proportions on the x axis. For visualization purposes,
                                  ## limit to the 0% spike-in, the 5 below the first that is significant,
                                  ## and all above the first that is significant
                                  limit.x.axis <- FALSE
                                  if(limit.x.axis) {
                                      orig.cmps <- cmps
                                      keep.indices <- c(1, max(1, min(sig.cmps$indx)-5):min(sig.cmps$indx), min(sig.cmps$indx):max(cmps$indx))
                                      keep.cmps <- orig.cmps[orig.cmps$indx %in% unique(keep.indices),]
                                      flag <- as.numeric(as.character(df.ds$measured)) %in% c(0, as.numeric(as.character(keep.cmps$group2)))
                                      df.ds <- df.ds[flag,]
                                  }
                                  
                                  g <- ggplot(data = df.ds, aes(x = measured, y = prediction))
                                  g <- g + geom_boxplot()
                                  ## Remove label = ..p.signif.. to plot the actual pvalues
                                  g <- g + stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, size = 5, vjust = 0.5)
                                  ## g <- g + geom_beeswarm()
                                  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                                  g <- g + xlab("Spike-in Percent") + ylab("Score")                                  
                                  ## g <- g + xlab("Spike-in Proportion") + ylab("Score")
                                  mt <- as.character(df.ds$mixture.type[1])
                                  title <- paste0(mt, " Admixtures")
                                  g <- g + ggtitle(title)
                                  g
                              })
                    ct <- as.character(df$cell.type[1])
                    meth <- as.character(df$repo_name[1])
                    subplot.names <- names(subplots)
                    names(subplot.names) <- subplot.names
                    l_ply(subplot.names,
                          .fun = function(sb.nm) {
                              mt <- sb.nm
                              prefix <- mt
                              title <- paste0(mt, " Admixtures")
                              g <- subplots[[sb.nm]]
                              g <- g + ggtitle(paste0(firstup(nm), "-Grained Sub-Challenge (Validation; ", title, "):\n", ct, " ", meth))
                              g <- g + theme(text = element_text(size = 20), title = element_text(size = 20))
                              file <- paste0("validation-lod-ct-", prefix, "-", ct, "-", meth, "-", nm)
                              file <- gsub(file, pattern = "\\.", replacement = "-")
                              file <- paste0(file, ".png")
                              png(file, width = 2 * 480)
                              print(g)
                              d <- dev.off()
                    })
                    
                    top <- paste0(firstup(nm), "-Grained Sub-Challenge (Validation):\n", ct, " ", meth)
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

