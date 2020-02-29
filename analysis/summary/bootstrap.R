library(ggplot2)

res <- read.table("results.csv", sep=",", header=TRUE, as.is=TRUE)
lb2 <- read.table("lb-r2.tsv", sep="\t", header=TRUE, as.is=TRUE)

res <- subset(res, ( is_latest == TRUE ) | ( grepl(repo_name, pattern = "baseline")))

res <- subset(res, round == 2)

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
                                                                         val <- cor(df.ct$prediction, df.ct$measured)
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

means.by.cell.type.method <-
    llply(bootstrapped.cors,
          .fun = function(df) {
              ## first, average over bootstrap
              ret <- ddply(df, .variables = c("method", "cell.type", "dataset.name"),
                           .fun = function(df) {
                               data.frame(cor = mean(df$cor, na.rm=TRUE))
                           })
              ## now, average over dataset
              ret2 <- ddply(ret, .variables = c("method", "cell.type"),
                            .fun = function(df) {
                                data.frame(cor = mean(df$cor, na.rm=TRUE))
                            })
              })

bootstrapped.scores <-
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
                                                                         val <- cor(df.ct$prediction, df.ct$measured)
                                                                         if(var(df.ct$prediction) == 0) { val <- 0 }
                                                                         data.frame(cor = val)
                                                                     })
                                                        if(FALSE) {
                                                        cell.types <- coarse.cell.types
                                                        if(sub.challenge == "fine") { cell.types <- fine.cell.types }
                                                        tmp <- ldply(cell.types, 
                                                                     .fun = function(ct) {
                                                                         if(!(ct %in% df.ds$cell.type)) {
                                                                             return(data.frame(cor = NA))
                                                                         }
                                                                         df.ct <- subset(df.ds, cell.type == ct)
                                                                         val <- cor(df.ct$prediction, df.ct$measured)
                                                                         if(var(df.ct$prediction) == 0) { val <- 0 }
                                                                         data.frame(cor = val)
                                                                     })
                                                        }
                                                        data.frame(score = mean(tmp$cor))
                                                        
                                                    })
                                          ret.method <- data.frame(score = mean(score$score))
                                          ret.method
                                      })
                            colnames(ret.all)[1] <- "method"
                            ret.all
                        })
          })

df <- bootstrapped.scores[["coarse"]]
means <- ddply(df, .variables = c("method"), .fun = function(tmp) data.frame(mean = mean(tmp$score, na.rm=TRUE)))
means <- means[order(means$mean),]
bootstrap.method.levels <- unique(means$method)
df$method <- factor(df$method, levels = bootstrap.method.levels)
g1 <- ggplot(data = df)
g1 <- g1 + geom_boxplot(aes(x = method, y = score))
g1 <- g1 + coord_flip()
g1 <- g1 + ggtitle("Coarse-Grained Sub-Challenge (Round 2)")
g1 <- g1 + ylab("Pearson Correlation\nAvg over Dataset & Cell Type")
g1 <- g1 + theme(text = element_text(size=20))
g1 <- g1 + xlab("")

df <- bootstrapped.scores[["fine"]]
means <- ddply(df, .variables = c("method"), .fun = function(tmp) data.frame(mean = mean(tmp$score, na.rm=TRUE)))
means <- means[order(means$mean),]
df$method <- factor(df$method, levels = unique(means$method))
g2 <- ggplot(data = df)
g2 <- g2 + geom_boxplot(aes(x = method, y = score))
g2 <- g2 + coord_flip()
g2 <- g2 + ggtitle("Fine-Grained Sub-Challenge (Round 2)")
g2 <- g2 + ylab("Pearson Correlation\nAvg over Dataset & Cell Type")
g2 <- g2 + theme(text = element_text(size=20))
g2 <- g2 + xlab("")

library(gridExtra)
library(cowplot)
## print(grid.arrange(g1,g2,nrow=1))
g <- plot_grid(g1, g2)
ggsave("round2-results.pdf", width = 14)

print(g)
ggsave("round2-results.png", width = 14)

pdf("round2-coarse-results.pdf")
print(g1)
d <- dev.off()

png("round2-coarse-results.png")
print(g1)
d <- dev.off()

tmp <- means.by.cell.type.method[["coarse"]]
library(dplyr)
library(tidyr)

tmp = tmp %>%
     group_by(method) %>%
     complete(cell.type = unique(tmp$cell.type))


means <- ddply(tmp, .variables = c("cell.type"),
               .fun = function(df) mean(df$cor, na.rm=TRUE))
means <- means[order(means$V1),]
tmp$cell.type <- factor(tmp$cell.type, levels = means$cell.type)

means <- ddply(tmp, .variables = c("method"),
               .fun = function(df) mean(df$cor, na.rm=TRUE))
means <- means[order(means$V1),]
## tmp$method <- factor(tmp$method, levels = means$method)

tmp$method <- factor(tmp$method, levels = bootstrap.method.levels)


g.coarse <- ggplot(data = tmp, 
                   aes(y = method, x = cell.type, fill = cor))
g.coarse <- g.coarse + geom_tile()
g.coarse <- g.coarse + theme(axis.text.x = element_text(angle = 90))
## g.coarse <- g.coarse + scale_fill_continuous("Pearson\ncorrelation", limits = c(-1,1))
g.coarse <- g.coarse + scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
g.coarse <- g.coarse + theme(text = element_text(size=20))
##g.coarse <- g.coarse + geom_point(data = tmp, aes(size="NA"), shape =NA, colour = "black")+
##    guides(size=guide_legend("", override.aes=list(shape=15, size = 10)))

png("round2-coarse-cell-type-results.png")
print(g.coarse)
d <- dev.off()
