suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))

suppressPackageStartupMessages(p_load("openxlsx"))
suppressPackageStartupMessages(p_load("reshape2"))


## Get the validation results ("all_predictions.csv")
synLogin()
synId <- "syn21715094"
## validation.results.file <- "validation-results.csv"
validation.results.file <- synGet(synId, downloadFile = TRUE)$path

res <- read.table(validation.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
## print(head(subset(res, method == "cibersort" & subchallenge == "fine" & cell.type == "fibroblasts" & dataset.name == "DS5")))

source("../utils.R")
val.metadata <- get.validation.metadata()
## flag <- !is.na(val.metadata$mixture.type) & (val.metadata$mixture.type == "BM")
## val.metadata$mixture.type[flag] <- "Biological"
## flag <- !is.na(val.metadata$mixture.type) & (val.metadata$mixture.type == "RM")
## val.metadata$mixture.type[flag] <- "Random"
    
res <- merge(res, val.metadata, all.x = TRUE, by.x = c("sample.id", "dataset.name"), by.y = c("id", "dataset"))

## res <- read.table("validation-results.csv", sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

## Limit to the final (of two) submissions for each group
## and the baseline methods (which were all submitted by Andrew L
## and hence only one of which is_latest)
## res <- subset(res, ( is_latest == TRUE ) | ( grepl(repo_name, pattern = "baseline")))
## res <- subset(res, !is.na(measured))
## Unlike for the leaderboard data, validation data is only composed of measured cell types.
## So, if something is NA, that means it was not present.
flag <- is.na(res$measured)
res[flag, "measured"] <- 0

exclude.purified <- TRUE
if(exclude.purified) {
    flag <- grepl(res$sample.id, pattern="BM") | grepl(res$sample.id, pattern="RM")
    res <- res[flag,]
    if(any(res$dataset.name == "DS5")) {
        stop("Was not expecting any admixtures (BM or RM) in dataset DS5\n")
    }
    ## Exclude neutrophils, which were not in the admixtures -- these just show up as having 0 measured.
    ## This creates problems. Not only do we need to add an epsilon, but methods actually score it well with
    ## low values, which will be uncorrelated since magnitude is not taken into account when during a correlation.
    ## Neither were naive.CD8.T.cells or memory.B.cells
    flag <- res$cell.type %in% c("neutrophils", "memory.B.cells", "naive.CD8.T.cells")
    res <- res[!flag, ]
    ## exclude <- c("memory.B.cells", "naive.CD8.T.cells")
    ## if(any(res$cell.type %in% exclude)) {
    ##    print(head(res[res$cell.type %in% exclude,]))
    ##    stop("Was not expecting memory.B.cells or naive.CD8.T.cells in admixtures\n")
    ## }
}

trans <-
    list("cibersort" = "CIBERSORT",
         "quantiseq" = "quanTIseq",
         "mcpcounter" = "MCP-counter",
         "epic" = "EPIC",
         "timer" = "TIMER",
         "xcell" = "xCell")
for(nm in names(trans)) {
    flag <- res$method == nm
    res[flag, "method"] <- trans[[nm]]
}

res$modelId <- res$method

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")

n.bootstraps <- 100
bootstraps <- 1:n.bootstraps
names(bootstraps) <- bootstraps

set.seed(1234)

bootstrapped.cors <-
    ddply(res,
          .variables = c("subchallenge"),
          .fun = function(df.sc) {
                    tmp <-
                        ldply(bootstraps,
                          .fun = function(i) {
                              ret <-
                                  ddply(df.sc, .variables = c("dataset.name", "mixture.type"),
                                        .fun = function(df.ds) {
                                            sample.ids <- unique(df.ds$sample.id)
                                            ## NB: use the same bootstrap samples across datasets and celltypes (nested below)
                                            boot.samples <- sample(sample.ids, size = length(sample.ids), replace = TRUE)
                                            eps.vec <- runif(length(boot.samples), min=10^-10, max=10^-9)
                                            ret.ds <-
                                                ddply(df.ds, .variables = c("modelId"),
                                                      .fun = function(df.method) {
                                                          ret.ct <-
                                                              ddply(df.method, .variables = c("cell.type"),
                                                                    .fun = function(df.ct) {
                                                                        rownames(df.ct) <- df.ct$sample.id
                                                                        ## This condition is necessary because a few datasets lack
                                                                        ## _measurements_ (and hence predictions) for some cell types in
                                                                        ## some (but not all) samples
                                                                        ev <- eps.vec[boot.samples %in% rownames(df.ct)]
                                                                        boot.samples <- boot.samples[boot.samples %in% rownames(df.ct)]
                                                                        if(length(boot.samples) < 2) { return(NULL) }
                                                                        boot.df <- df.ct[boot.samples, c("prediction", "measured")]
                                                                        cor.p <- 0
                                                                        cor.s <- 0
                                                                        cor.p.pval <- NA
                                                                        cor.s.pval <- NA
                                                                        ## In some cases all measured values are zero, if so just add
                                                                        ## a little noise
                                                                        if(!all(boot.df$measured == boot.df$measured[1])) { ev <- 0 }
                                                                        if(var(boot.df$prediction) != 0) {
                                                                            cor.p.out <- cor.test(boot.df$prediction, boot.df$measured + ev,
                                                                                                  method = "pearson")
                                                                            cor.p <- as.numeric(cor.p.out$estimate)
                                                                            cor.p.pval <- as.numeric(cor.p.out$p.value)
                                                                            if(is.na(cor.p)) {
                                                                                print(boot.df)
                                                                                stop("stop")
                                                                            }
                                                                            cor.s.out <- cor.test(boot.df$prediction, boot.df$measured + ev,
                                                                                                  method = "spearman")
                                                                            cor.s <- as.numeric(cor.s.out$estimate)
                                                                            cor.s.pval <- as.numeric(cor.s.out$p.value)
                                                                            if(is.na(cor.s)) {
                                                                                print(boot.df)
                                                                                stop("stop")
                                                                            }
                                                                        }
                                                                        data.frame(cor.p = cor.p, cor.p.pval = cor.p.pval,
                                                                                   cor.s = cor.s, cor.s.pval = cor.s.pval)
                                                                    }) # function(df.ct)
                                                      }) # function(df.method)
                                        }) # function(df.ds)
                              ret
                          }) # function(i)
                    colnames(tmp)[1] <- "bootstrap"
                    tmp
          })

write.table(file = "validation-bootstraps.tsv", bootstrapped.cors, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## bootstrapped.cors <- read.table("bootstraps.tsv", sep = "\t", header = TRUE)

bootstrapped.scores <-
    ddply(bootstrapped.cors,
          .variables = c("subchallenge", "bootstrap", "modelId"),
          .fun = function(df) {
              res.ds <-
                  ddply(df, .variables = c("dataset.name"),
                        .fun = function(df.ds) {
                            data.frame(ds.score = mean(df$cor.p))
                        })
              data.frame(score = mean(res.ds$ds.score))
          })

bootstrapped.scores.sv <- bootstrapped.scores

plot.bootstraps <- function(df) {
    means <- ddply(df, .variables = c("modelId"),
                   .fun = function(tmp) data.frame(mean = mean(tmp$score, na.rm=TRUE)))
    means <- means[order(means$mean),]
    means <- means[!duplicated(means$modelId, fromLast = TRUE),]
    df$modelId <- factor(df$modelId, levels = unique(means$modelId))
    g1 <- ggplot(data = df)
    g1 <- g1 + geom_boxplot(aes(x = modelId, y = score))
    int.data <- data.frame(mean = mean(means$mean))
    g1 <- g1 + geom_hline(data = int.data, aes(yintercept = mean), linetype = "dashed")
    g1 <- g1 + coord_flip()
    g1 <- g1 + ylab("Pearson Correlation\nAvg over Dataset & Cell Type")
    g1 <- g1 + theme(text = element_text(size=20), 
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     panel.spacing = unit(1, "lines"))
    g1 <- g1 + xlab("Method")
    g1
}

mean.bootstrapped.scores <-
    ddply(bootstrapped.scores,
          .variables = c("subchallenge", "modelId"),
          .fun = function(df) {
              data.frame(mean.score = mean(df$score))
          })

means.by.cell.type.method <-
    ddply(bootstrapped.cors,
          .variables = c("subchallenge", "modelId", "cell.type"),
          .fun = function(df) {
              ## first, average over bootstrap
              ret <- ddply(df, .variables = c("dataset.name"),
                           .fun = function(df) {
                               data.frame(cor.p = mean(df$cor.p, na.rm=TRUE), cor.p.pval = mean(df$cor.p.pval, na.rm=TRUE),
                                          cor.s = mean(df$cor.s, na.rm=TRUE), cor.s.pval = mean(df$cor.s.pval, na.rm=TRUE))
                           })
              ## now, average over dataset
              data.frame(cor.p = mean(ret$cor.p), cor.p.pval = mean(ret$cor.p.pval),
                         cor.s = mean(ret$cor.s), cor.s.pval = mean(ret$cor.s.pval))
          })

means.by.cell.type.method.dataset <-
    ddply(bootstrapped.cors,
          .variables = c("subchallenge", "modelId", "cell.type", "dataset.name", "mixture.type"),
          .fun = function(df) {
              ## average over bootstrap
              data.frame(cor.p = mean(df$cor.p), cor.p.pval = mean(df$cor.p.pval),
                         cor.s = mean(df$cor.s), cor.s.pval = mean(df$cor.s.pval))
          })

means.by.cell.type.method.mixture.type <-
    ddply(bootstrapped.cors,
          .variables = c("subchallenge", "modelId", "cell.type", "mixture.type"),
          .fun = function(df) {
              ## first, average over bootstrap
              ret <- ddply(df, .variables = c("dataset.name"),
                           .fun = function(df) {
                               data.frame(cor.p = mean(df$cor.p, na.rm=TRUE), cor.p.pval = mean(df$cor.p.pval, na.rm=TRUE),
                                          cor.s = mean(df$cor.s, na.rm=TRUE), cor.s.pval = mean(df$cor.s.pval, na.rm=TRUE))
                           })
              ## now, average over dataset
              data.frame(cor.p = mean(ret$cor.p), cor.p.pval = mean(ret$cor.p.pval),
                         cor.s = mean(ret$cor.s), cor.s.pval = mean(ret$cor.s.pval))
          })



plot.cell.type.correlation.heatmap.old <- function(df) {
    means <- ddply(df, .variables = c("cell.type"),
                   .fun = function(tmp) mean(tmp$cor, na.rm=TRUE))
    means <- means[order(means$V1),]
    df$cell.type <- factor(df$cell.type, levels = means$cell.type)
    
    means <- ddply(df, .variables = c("modelId"),
                   .fun = function(tmp) mean(tmp$cor, na.rm=TRUE))
    means <- means[order(means$V1),]
    df$modelId <- factor(df$modelId, levels = means$modelId)
    df$cor.label <- formatC(df$cor, format="f", digits=2)
    g <- ggplot(data = df, aes(y = modelId, x = cell.type, fill = cor))
    g <- g + geom_tile()
    g <- g + geom_text(aes(label = cor.label))
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   title = element_text(size = 8))
    g <- g + ylab("Method") + xlab("")
    ## g <- g + scale_fill_continuous("Pearson\ncorrelation", limits = c(-1,1))
    ## g <- g + scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1),
    ##                               low = "red", high = "blue", mid = "white", na.value = "black")
    g <- g + scale_fill_gradient2("Correlation", limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
    ## g <- g + theme(text = element_text(size=20))
    g
}

df <- subset(bootstrapped.scores, subchallenge == "coarse")
g1 <- plot.bootstraps(df)
g1 <- g1 + ggtitle("Coarse-Grained Sub-Challenge\n(Validation)")
## pdf("validation-coarse-bootstrap.pdf", width = 14)
png("validation-coarse-bootstrap.png", width = 480 * 1)
print(g1)
d <- dev.off()

df <- subset(bootstrapped.scores, subchallenge == "fine")
g2 <- plot.bootstraps(df)
g2 <- g2 + ggtitle("Fine-Grained Sub-Challenge\n(Validation)")
## pdf("validation-fine-bootstrap.pdf", width = 14)
png("validation-fine-bootstrap.png", width = 480 * 1)
print(g2)
d <- dev.off()

png("validation-coarse-cell-type.png", width = 1 * 480)
g3 <- plot.cell.type.correlation.heatmap(subset(means.by.cell.type.method, subchallenge == "coarse"),
                                         show.corr.text = TRUE, id.var = "modelId")
g3 <- g3 + ggtitle("Coarse-Grained Sub-Challenge (Validation)")
print(g3)
d <- dev.off()

png("validation-fine-cell-type.png", width = 1 * 480)
g4 <- plot.cell.type.correlation.heatmap(subset(means.by.cell.type.method, subchallenge == "fine"),
                                         show.corr.text = TRUE, id.var = "modelId")
g4 <- g4 + ggtitle("Fine-Grained Sub-Challenge (Validation)")
print(g4)
d <- dev.off()


cor.types <- list("Pearson" = "Pearson", "Spearman" = "Spearman")
ds.plts <-
    llply(cor.types,
          .fun = function(cor.type) {
              dlply(means.by.cell.type.method.dataset,
                    .variables = c("subchallenge", "dataset.name"),
                    .fun = function(df) {
                        cor.var <- "cor.p"
                        pval.var <- "cor.p.pval"
                        if(cor.type == "Spearman") { cor.var <- "cor.s"; pval.var <- "cor.s.pval" }
                        g <- plot.cell.type.correlation.heatmap(df, show.corr.text = TRUE, id.var = "modelId",
                                                                cor.var = cor.var, cor.type.label = paste0(cor.type, "\nCorrelation"),
                                                                pval.var = pval.var)
                        sc <- as.character(df$subchallenge[1])
                        ds <- as.character(df$dataset.name[1])
                        mt <- as.character(df$mixture.type[1])
                        g <- g + ggtitle(paste0(firstup(sc), "-Grained Sub-Challenge\n(Validation; ", ds, "; ", mt, ")"))
                        g
                    })
          })

mt.plts <-
    llply(cor.types,
          .fun = function(cor.type) {
              dlply(means.by.cell.type.method.mixture.type,
                    .variables = c("subchallenge", "mixture.type"),
                    .fun = function(df) {
                        cor.var <- "cor.p"
                        pval.var <- "cor.p.pval"
                        if(cor.type == "Spearman") { cor.var <- "cor.s"; pval.var <- "cor.s.pval" }
                        g <- plot.cell.type.correlation.heatmap(df, show.corr.text = TRUE, id.var = "modelId",
                                                                cor.var = cor.var, cor.type.label = paste0(cor.type, "\nCorrelation"),
                                                                pval.var = pval.var)
                        sc <- as.character(df$subchallenge[1])
                        mt <- as.character(df$mixture.type[1])
                        g <- g + ggtitle(paste0(firstup(sc), "-Grained Sub-Challenge\n(Validation; ", mt, ")"))
                        g
                    })
          })

for(cor.type in c("pearson", "spearman")) {
    for(sc in c("coarse", "fine")) {
        g1 <- mt.plts[[firstup(cor.type)]][[paste0(sc, ".Random")]]
        g1 <- g1 + ggtitle("Random Admixtures")
        g2 <- mt.plts[[firstup(cor.type)]][[paste0(sc, ".Biological")]]
        g2 <- g2 + ggtitle("Biological Admixtures")
        title <- paste0(firstup(sc), "-Grained Sub-Challenge (Validation)")
        png(paste0("validation-", sc, "-bio-and-rand-", cor.type, ".png"), width = 2 * 480)
        g <- grid.arrange(g1, g2, ncol = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
        grid.draw(g)
        d <- dev.off()
    }
}


## Collect the ground truth data
##gt <- unique(res[, c("sample.id", "subchallenge", "dataset.name", "mixture.type", "tumor.type", "cell.type", "measured")])
gt <- get.validation.ground.truth()

plts <- 
    dlply(gt,
          .variables = c("mixture.type"),
          .fun = function(df) {
              mt <- as.character(df$mixture.type[1])
              ds <- as.character(df$dataset.name[1])
              df <- df[, c("sample.id", "cell.type", "measured")]
              mat <- acast(df, cell.type ~ sample.id)
              g <- plot.admixtures(mat)
              ## g <- g + ggtitle(paste0(firstup(sc), "-Grained Sub-Challenge\n(Validation; ", ds, "; ", mt, ")"))
              ## g <- g + ggtitle(paste0("Validation Ground Truth: ", ds, ", ", mt))
              g <- g + ggtitle(paste0("Validation Ground Truth:\n", mt, " Admxitures"))
              g
          })

d_ply(res,
      .variables = c("cell.type", "subchallenge", "modelId"),
      .fun = function(df) {
          ct <- as.character(df$cell.type[1])
          sc <- as.character(df$subchallenge[1])
          meth <- as.character(df$modelId[1])
          mt <- "Biological"
          sub <- subset(df, mixture.type == mt)

          g1 <- ggplot(sub, aes(x = measured, y = prediction))
          g1 <- g1 + geom_point()
          g1 <- g1 + facet_wrap(~ dataset.name, scales = "free")
          g1 <- g1 + geom_smooth(method = "lm")
          g1 <- g1 + ggtitle("Biological Admixtures")
          g1 <- g1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
          
          mt <- "Random"
          sub <- subset(df, mixture.type == mt)          
          
          g2 <- ggplot(sub, aes(x = measured, y = prediction))
          g2 <- g2 + geom_point()
          g2 <- g2 + facet_wrap(~ dataset.name, scales = "free")
          g2 <- g2 + geom_smooth(method = "lm")
          g2 <- g2 + ggtitle("Random Admixtures")
          g2 <- g2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
          
          title <- paste0(meth, " ", ct, " (", firstup(sc), "-Grained Sub-Challenge)")
          g <- grid.arrange(g2, g1, ncol = 2, top = textGrob(title, gp = gpar(fontsize = 25)))

          png(paste0("validation-admixtures-", "method-", meth, "-", ct, "-", sc, ".png"), width = 2 * 480)
          grid.draw(g)
          d <- dev.off()
      })



png(paste0("validation-admixtures", ".png"), width = 2 * 480)
g1 <- plts[["Random"]]
g1 <- g1 + ggtitle("Random Admixtures")
g2 <- plts[["Biological"]]
g2 <- g2 + ggtitle("Biological Admixtures")
title <- "Validation Admixture Proportions"
g <- grid.arrange(g1, g2, ncol = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
grid.draw(g)
d <- dev.off()
