library(ggplot2)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

res <- read.table("results.csv", sep=",", header=TRUE, as.is=TRUE)

res <- subset(res, ( is_latest == TRUE ) | ( grepl(repo_name, pattern = "baseline")))

res <- subset(res, round == 2)


## Score using Aurelien's approach -- i.e., weighting by the std dev of correlations
## (over methods)

res$modelId <- paste0(as.character(res$repo_name), "-", as.character(res$submitterId))

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")

library(plyr)
library(dplyr)

tbls <-
    llply(sub.challenges,
          .fun = function(sub.challenge) {
              tmp <- subset(res, subchallenge == sub.challenge)
              tmp <- subset(tmp, !is.na(measured))
              tmp
          })

tbls.by.cell <-
    llply(sub.challenges,
          .fun = function(sub.challenge) {
              tmp <- subset(res, subchallenge == sub.challenge)
              tmp <- subset(tmp, !is.na(measured))
              cells <- unique(tmp$cell.type)
              names(cells) <- cells
              llply(cells,
                    .fun = function(ct) {
                        ret <- subset(tmp, cell.type == ct)
                        ret$id <- paste0(ret$dataset.name, "-", ret$sample.id)
                        ret
                    })
          })
    
coarse.cell.types <- unique(tbls[["coarse"]]$cell.type)
fine.cell.types <- unique(tbls[["fine"]]$cell.type)

n.bootstraps <- 100

bootstraps <-
    llply(sub.challenges,
          .fun = function(sub.challenge) {
              tmp <- subset(res, subchallenge == sub.challenge)
              tmp <- subset(tmp, !is.na(measured))
              datasets <-
                  dlply(tmp, .variables = c("dataset.name"),
                        .fun = function(df) {
                            unique(df$sample.id)
                        })
              dataset.names <- names(datasets)
              names(dataset.names) <- dataset.names
              llply(1:n.bootstraps,
                    .fun = function(i) {
                        ret <- ldply(datasets,
                                     .fun = function(ds) {
                                         data.frame(sample.id = sample(ds, size = length(ds), replace = TRUE))
                                         
                                     })
                        colnames(ret)[1] <- "dataset.name"
                        ret$id <- paste0(ret$dataset.name, "-", ret$sample.id)                        
                        ret
                    })
          })

bootstrapped.cors <-
    llply(sub.challenges,
          .fun = function(sub.challenge) {
              tmp <- subset(res, subchallenge == sub.challenge)
              tmp <- subset(tmp, !is.na(measured))
##              tmp$modelId <- paste0(tmp$repo_name, "-", tmp$submitterId)
              ##              methods <- unique(tmp$repo_name)
              methods <- unique(tmp$modelId)
              names(methods) <- methods
              indices <- 1:n.bootstraps
              names(indices) <- indices
              ret.i <-
                  ldply(indices,
                        .fun = function(i) {
                            ret.all <-
                                ldply(methods,
                                      .fun = function(method) {
                                          ret.method <-
                                              ldply(tbls.by.cell[[sub.challenge]],
                                                    .fun = function(df) {
                                                        ##                                                        df <- subset(df, repo_name == method)
                                                        df <- subset(df, modelId == method)
                                                        rownames(df) <- df$id
                                                        sample.ids <- bootstraps[[sub.challenge]][[i]]$id
                                                        sample.ids <- sample.ids[sample.ids %in% df$id]
                                                        df <- df[sample.ids,]
                                                        df
                                                    })
                                          colnames(ret.method)[1] <- method
                                          score <-
                                              ddply(ret.method, .variables = c("dataset.name"),
                                                    .fun = function(df.ds) {
                                                        tmp <- ddply(df.ds, .variables = c("cell.type"),
                                                                     .fun = function(df.ct) {
                                                                         val <- cor(df.ct$prediction, df.ct$measured, method = "pearson")
                                                                         if(var(df.ct$prediction) == 0) { val <- 0 }
                                                                         data.frame(cor = val)
                                                                     })
                                                        colnames(tmp)[1] <- "cell.type"
                                                        tmp
                                                    })
                                          score
                                      })
                            colnames(ret.all)[1] <- "method"
                            ret.all
                        })
              colnames(ret.i)[1] <- "boot.i"
              ret.i
          })

scores.by.cell.type.method.boot <-
    llply(bootstrapped.cors,
          .fun = function(tmp.df) {
              ## Define std dev over correlation within each bootstrap
              ret1 <- ddply(tmp.df, .variables = c("boot.i", "cell.type", "dataset.name"),
                            .fun = function(df) {
                                std.dev <- sd(df$cor, na.rm=TRUE)
                                data.frame(cor.sd = std.dev)
                            })
              ret2 <- merge(ret1, tmp.df)

              ## Define score as correlation weighted by standard deviation
              ret2$score <- ret2$cor.sd * ret2$cor
              
              ## Sum score over dataset (within each bootstrap)
              ret3 <- ddply(ret2, .variables = c("boot.i", "method", "cell.type"),
                            .fun = function(df) {
                                data.frame(score = sum(df$score, na.rm=TRUE))
                            })
          })

mean.correlations.by.cell.type.method.boot <-
    llply(bootstrapped.cors,
          .fun = function(tmp.df) {
              ## Take mean correlation over dataset (within each bootstrap)
              ret3 <- ddply(tmp.df, .variables = c("boot.i", "method", "cell.type"),
                            .fun = function(df) {
                                data.frame(score = mean(df$cor, na.rm=TRUE))
                            })
          })


## Sum score over cell type (within each bootstrap)
scores.by.method.boot <-
    llply(scores.by.cell.type.method.boot,
          .fun = function(tmp) {
              ret <- ddply(tmp, .variables = c("boot.i", "method"),
                           .fun = function(df) {
                               data.frame(score = sum(df$score, na.rm=TRUE))
                           })
          })

## Mean correlation over cell type (within each bootstrap)
mean.correlations.by.method.boot <-
    llply(mean.correlations.by.cell.type.method.boot,
          .fun = function(tmp) {
              ret <- ddply(tmp, .variables = c("boot.i", "method"),
                           .fun = function(df) {
                               data.frame(score = mean(df$score, na.rm=TRUE))
                           })
          })

## Define mean score (over bootstraps) by cell type / method
mean.correlations.by.cell.type.method <-
    llply(mean.correlations.by.cell.type.method.boot,
          .fun = function(tmp) {
              ret <- ddply(tmp, .variables = c("cell.type", "method"),
                           .fun = function(df) {
                               data.frame(score = mean(df$score, na.rm=TRUE))
                           })
          })
              
library(gridExtra)
library(cowplot)


approach <- "sd-weighted-sum"

for(sc in sub.challenges) {

    df <- scores.by.method.boot[[sc]]
    means <- ddply(df, .variables = c("method"), .fun = function(tmp) data.frame(mean = mean(tmp$score, na.rm=TRUE)))
    means <- means[order(means$mean),]
    bootstrap.method.levels <- unique(means$method)
    df$method <- factor(df$method, levels = bootstrap.method.levels)
    g1 <- ggplot(data = df)
    g1 <- g1 + geom_boxplot(aes(x = method, y = score))
    g1 <- g1 + coord_flip()
    g1 <- g1 + ggtitle(paste0(firstup(sc), "-Grained (Round 2)"))
    g1 <- g1 + ylab("Std Dev-Weighted\nPearson Correlation")
    g1 <- g1 + theme(text = element_text(size=20))
    g1 <- g1 + xlab("")
    
    ##g <- plot_grid(g1, g2)
    ##ggsave(paste0("round2-results-", approach, ".pdf"), width = 14)
    
    ##print(g)
    ##ggsave(paste0("round2-results-", approach, ".png"), width = 14)
    
    pdf(paste0("round2-", sc, "-results-", approach, ".pdf"))
    print(g1)
    d <- dev.off()
    
    png(paste0("round2-", sc, "-results-", approach, ".png"))
    print(g1)
    d <- dev.off()
    
    tmp <- mean.correlations.by.cell.type.method[[sc]]
    library(dplyr)
    library(tidyr)
    
    tmp = tmp %>%
        group_by(method) %>%
        complete(cell.type = unique(tmp$cell.type))
    
    
    means <- ddply(tmp, .variables = c("cell.type"),
                   .fun = function(df) mean(df$score, na.rm=TRUE))
    means <- means[order(means$V1),]
    sds <- ddply(tmp, .variables = c("cell.type"),
                 .fun = function(df) sd(df$score, na.rm=TRUE))
    sds <- sds[order(sds$V1),]
    ## tmp$cell.type <- factor(tmp$cell.type, levels = means$cell.type)
    tmp$cell.type <- factor(tmp$cell.type, levels = sds$cell.type)
    
    ## means <- ddply(tmp, .variables = c("method"),
    ##                .fun = function(df) mean(df$score, na.rm=TRUE))
    ##  means <- means[order(means$V1),]
    ## tmp$method <- factor(tmp$method, levels = means$method)
    
    tmp$method <- factor(tmp$method, levels = bootstrap.method.levels)
    
    
    g.plt <- ggplot(data = tmp, 
                    aes(y = method, x = cell.type, fill = score))
    g.plt <- g.plt + xlab("Cell Type\nOrdered by Std Dev")
    g.plt <- g.plt + ylab("Method\nOrdered by Std Dev-Weighted Correlation")
    g.plt <- g.plt + ggtitle(paste0(firstup(sc), "-Grained (Round 2)"))
    g.plt <- g.plt + geom_tile()
    g.plt <- g.plt + theme(axis.text.x = element_text(angle = 90))
    ## g.plt <- g.plt + scale_fill_continuous("Pearson\ncorrelation", limits = c(-1,1))
    g.plt <- g.plt + scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
    g.plt <- g.plt + theme(text = element_text(size=20))
    ##g.plt <- g.plt + geom_point(data = tmp, aes(size="NA"), shape =NA, colour = "black")+
    ##    guides(size=guide_legend("", override.aes=list(shape=15, size = 10)))
    
    png(paste0("round2-", sc, "-cell-type-results-", approach, ".png"))
    print(g.plt)
    d <- dev.off()
}

## Repeat above using weighted correlation
approach <- "weighted-correlation"

for(sc in sub.challenges) {
    df <- mean.correlations.by.method.boot[[sc]]
    means <- ddply(df, .variables = c("method"), .fun = function(tmp) data.frame(mean = mean(tmp$score, na.rm=TRUE)))
    means <- means[order(means$mean),]
    bootstrap.method.levels <- unique(means$method)
    df$method <- factor(df$method, levels = bootstrap.method.levels)
    g1 <- ggplot(data = df)
    g1 <- g1 + geom_boxplot(aes(x = method, y = score))
    g1 <- g1 + coord_flip()
    g1 <- g1 + ggtitle(paste0(firstup(sc), "-Grained (Round 2)"))
    g1 <- g1 + ylab("Pearson Correlation\nAvg over Dataset then Cell Type")
    g1 <- g1 + theme(text = element_text(size=20))
    g1 <- g1 + xlab("")
    
    ##g <- plot_grid(g1, g2)
    ##ggsave(paste0("round2-results-", approach, ".pdf"), width = 14)
    
    ##print(g)
    ##ggsave(paste0("round2-results-", approach, ".png"), width = 14)
    
    pdf(paste0("round2-", sc, "-results-", approach, ".pdf"))
    print(g1)
    d <- dev.off()
    
    png(paste0("round2-", sc, "-results-", approach, ".png"))
    print(g1)
    d <- dev.off()
    
    tmp <- mean.correlations.by.cell.type.method[[sc]]
    library(dplyr)
    library(tidyr)
    
    tmp = tmp %>%
        group_by(method) %>%
        complete(cell.type = unique(tmp$cell.type))
    
    means <- ddply(tmp, .variables = c("cell.type"),
                   .fun = function(df) mean(df$score, na.rm=TRUE))
    means <- means[order(means$V1),]
    tmp$cell.type <- factor(tmp$cell.type, levels = means$cell.type)
    
    sds <- ddply(tmp, .variables = c("cell.type"),
                 .fun = function(df) sd(df$score, na.rm=TRUE))
    sds <- sds[order(sds$V1),]
    tmp$cell.type <- factor(tmp$cell.type, levels = sds$cell.type)
    
    ##means <- ddply(tmp, .variables = c("method"),
    ##               .fun = function(df) mean(df$score, na.rm=TRUE))
    ##means <- means[order(means$V1),]
    ## tmp$method <- factor(tmp$method, levels = means$method)
    
    tmp$method <- factor(tmp$method, levels = bootstrap.method.levels)
    
    g.plt <- ggplot(data = tmp, 
                    aes(y = method, x = cell.type, fill = score))
    g.plt <- g.plt + xlab("Cell Type\nOrdered by Std Dev")
    g.plt <- g.plt + ylab("Method\nOrdered by Avg Correlation")
    g.plt <- g.plt + ggtitle(paste0(firstup(sc), "-Grained (Round 2)"))    
    g.plt <- g.plt + geom_tile()
    g.plt <- g.plt + theme(axis.text.x = element_text(angle = 90))
    ## g.plt <- g.plt + scale_fill_continuous("Pearson\ncorrelation", limits = c(-1,1))
    g.plt <- g.plt + scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
    g.plt <- g.plt + theme(text = element_text(size=20))
    ##g.plt <- g.plt + geom_point(data = tmp, aes(size="NA"), shape =NA, colour = "black")+
    ##    guides(size=guide_legend("", override.aes=list(shape=15, size = 10)))
    
    png(paste0("round2-", sc, "-cell-type-results-", approach, ".png"))
    print(g.plt)
    d <- dev.off()
}
