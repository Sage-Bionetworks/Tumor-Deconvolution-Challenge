suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))


## Get the validation results ("all_predictions.csv")
synLogin()
synId <- "syn21715094"
## validation.results.file <- "validation-results.csv"
validation.results.file <- synGet(synId, downloadFile = TRUE)$path

res <- read.table(validation.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
## print(head(subset(res, method == "cibersort" & subchallenge == "fine" & cell.type == "fibroblasts" & dataset.name == "DS5")))

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
                                  ddply(df.sc, .variables = c("dataset.name"),
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
                                                                        val <- 0
                                                                        ## In some cases all measured values are zero, if so just add
                                                                        ## a little noise
                                                                        if(!all(boot.df$measured == boot.df$measured[1])) { ev <- 0 }
                                                                        if(var(boot.df$prediction) != 0) {
                                                                            val <- cor(boot.df$prediction, boot.df$measured + ev)
                                                                            if(is.na(val)) {
                                                                                print(boot.df)
                                                                                stop("stop")
                                                                            }
                                                                        }
                                                                        data.frame(cor = val)
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
                            data.frame(ds.score = mean(df$cor))
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
                               data.frame(cor = mean(df$cor, na.rm=TRUE))
                           })
              ## now, average over dataset
              data.frame(cor = mean(ret$cor))
          })

plot.cell.type.correlation.heatmap <- function(df, show.corr.text = FALSE, id.var = "modelId") {

    df <- df[, c(id.var, "cell.type", "cor")]
    df[, id.var] <- as.character(df[, id.var])
    df$cell.type <- as.character(df$cell.type)
    
    cell.type.means <- ddply(df, .variables = c("cell.type"),
                             .fun = function(tmp) {
                                 ret <- data.frame(id = "mean", cor = mean(tmp$cor, na.rm=TRUE))
                                 colnames(ret)[1] <- id.var
                                 ret
                             })
    cell.type.means <- cell.type.means[order(cell.type.means$cor),]
    
    method.means <- ddply(df, .variables = c(id.var),
                   .fun = function(tmp) data.frame(cell.type = "mean", cor = mean(tmp$cor, na.rm=TRUE)))
    method.means <- method.means[order(method.means$cor),]


    cell.type.levels <- c(cell.type.means$cell.type[cell.type.means$cell.type != "mean"], "mean")
    id.levels <- c("mean", method.means[method.means[, id.var] != "mean", id.var])

    df <- rbind(df, cell.type.means[, c(id.var, "cell.type", "cor")])
    df <- rbind(df, method.means[, c(id.var, "cell.type", "cor")])    
    df$cell.type <- factor(df$cell.type, levels = cell.type.levels)
    df[, id.var] <- factor(df[, id.var], levels = id.levels)

    df$cor.label <- formatC(df$cor, format="f", digits=2)
    g <- ggplot(data = df, aes_string(y = id.var, x = "cell.type", fill = "cor"))
    g <- g + geom_tile()
    if(show.corr.text) {
        g <- g + geom_text(aes(label = cor.label))
    }
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   text = element_text(size=15))
    g <- g + ylab("Method") + xlab("")
    ## g <- g + scale_fill_continuous("Pearson\ncorrelation", limits = c(-1,1))
    ## g <- g + scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1),
    ##                               low = "red", high = "blue", mid = "white", na.value = "black")
    g <- g + scale_fill_gradient2("Correlation", limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
    ## g <- g + theme(text = element_text(size=20))
    g
}

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


