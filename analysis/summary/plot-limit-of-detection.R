
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

synLogin()
synId <- "syn21739521"
in.silico.results.file <- synGet(synId, downloadFile = TRUE)$path
res <- read.table(in.silico.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

synId <- "syn21752552"
in.silico.coarse.admixture.file <- synGet(synId, downloadFile = TRUE)$path
in.silico.coarse.admixtures <- read.table(in.silico.coarse.admixture.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

synId <- "syn21752551"
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

spike.in.res <- list()
spike.in.res[["coarse"]] <- subset(res, subchallenge == "coarse")
tmp <- na.omit(unique(in.silico.coarse.admixtures[, c("dataset.name", "sample.id", "spike.in.pop", "spike.in.prop")]))
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
              df$measured <- sum(df$measured)
              df[1,,drop=F]
          })
spike.in.res[["coarse"]] <- merge(spike.in.res[["coarse"]], tmp, by.x = c("dataset.name", "sample.id", "cell.type"),
                                  by.y = c("dataset.name", "sample.id", "spike.in.pop"))
                                  

spike.in.res[["fine"]] <- subset(res, subchallenge == "fine")
tmp <- na.omit(unique(in.silico.fine.admixtures[, c("dataset.name", "sample.id", "spike.in.pop", "spike.in.prop")]))
spike.in.res[["fine"]] <- merge(spike.in.res[["fine"]], tmp, by.x = c("dataset.name", "sample.id", "cell.type"),
                                  by.y = c("dataset.name", "sample.id", "spike.in.pop"))

res.matrices <-
    llply(spike.in.res,
          .fun = function(res) {
              ddply(res,
                    .variables = c("subchallenge", "cell.type", "dataset.name", "repo_name"),
                    .fun = function(df) {
                        df <- df[order(df$measured),]
                        sci <- formatC(as.numeric(as.character(df$measured)), format="e", digits=2)
                        df$measured <- factor(sci, levels = unique(sci))
                        
                        cmps <- as.data.frame(compare_means(prediction ~ measured,  data = df))
                        cmps <- cmps[cmps$group1 == "0.00e+00",]
                        ret <- cmps[, c("group1", "group2", "p")]
                        min.diff.prop <- 1
                        if(any(ret$p < 0.01)) {
                            flag <- ret$p < 0.01
                            min.diff.prop <- min(as.numeric(ret[flag,"group2"]))
                        }
                        data.frame(min.diff.prop = min.diff.prop)
                    })
          })

spike.ins <- c(1/(2^seq(from=2,to=12,by=1)), 0)
names(spike.ins) <- spike.ins

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

nms <- names(res.matrices)
names(nms) <- nms
plts <-
    llply(nms,
          .fun = function(nm) {
              res <- res.matrices[[nm]]
              subplots <-
                  dlply(res, .variables = c("dataset.name"),
                        .fun = function(res.ds) {

                            res.ds$label <- formatC(as.numeric(as.character(res.ds$min.diff.prop)), format="e", digits=1)
                            g <- ggplot(data = res.ds, aes(y = repo_name, x = cell.type, fill = log2(min.diff.prop)))
                            g <- g + geom_tile()
                            g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                           title = element_text(size = 8))
                            g <- g + ylab("Method") + xlab("")
                            g <- g + scale_fill_gradient2("LoD", limits = c(min(unname(log2(spike.ins[spike.ins>0]))),log2(1)), low = "red", high = "blue", mid = "white", na.value = "black")
                            g <- g + geom_text(aes(label = label))
                            title <- "blank"
                            if(res.ds[1,"dataset.name"] %in% c("A", "G")) {
                                title <- "Random Admixtures"
                            }
                            if(res.ds[1,"dataset.name"] %in% c("B", "H")) {
                                title <- "Biological Admixtures"
                            }
                            g <- g + ggtitle(title)
                        })
              do.call("grid.arrange", c(subplots, nrow = 2, top = paste0(firstup(nm), "-Grained Sub-Challenge (Validation)")))
          })

l_ply(nms,
      .fun = function(nm) {
          png(paste0("validation-lod-summary-", nm, ".png"), width = 2 * 480)
          grid.draw(plts[[nm]])
          d <- dev.off()
      })

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
                              .variables = c("dataset.name"),
                              .fun = function(df.ds) {

                                  df.ds <- df.ds[order(df.ds$measured),]
                                  sci <- formatC(as.numeric(as.character(df.ds$measured)), format="e", digits=2)
                                  df.ds$measured <- factor(sci, levels = unique(sci))
                                  cmps <- as.data.frame(compare_means(prediction ~ measured,  data = df.ds))
                                  cmps <- cmps[cmps$group1 == "0.00e+00",]
                                  my_comparisons <-
                                      llply(1:nrow(cmps),
                                            .fun = function(i) c(cmps[i,2], cmps[i,3]))
                                  
                                  
                                  g <- ggplot(data = df.ds, aes(x = measured, y = prediction))
                                  g <- g + geom_boxplot()
                                  g <- g + stat_compare_means(comparisons = my_comparisons)
                                  g <- g + geom_beeswarm()
                                  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                                  g <- g + xlab("Spike-in Proportion") + ylab("Score")
                                  title <- "blank"
                                  if(df.ds[1,"dataset.name"] %in% c("A", "G")) {
                                      title <- "Random Admixtures"
                                  }
                                  if(df.ds[1,"dataset.name"] %in% c("B", "H")) {
                                      title <- "Biological Admixtures"
                                  }
                                  g <- g + ggtitle(title)
                                  g
                              })
                    ct <- as.character(df$cell.type[1])
                    meth <- as.character(df$repo_name[1])
                    top <- paste0(firstup(nm), "-Grained Sub-Challenge (Validation): ", ct, " ", meth)
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

stop("stop")

sc <- "coarse"
## tmp <- subset(res, (measured %in% spike.ins) & (subchallenge == sc) & (dataset.name %in% c("G", "H")))
tmp <- subset(res, abs(measured - 0.25) < 10^-4 & (subchallenge == sc) & (dataset.name %in% c("G", "H")))
foo <- subset(tmp, measured > 0)
## Do this via subtraction / eps
tmp <- subset(tmp, cell.type %in% unique(foo$cell.type))
##flag <- (tmp$measured %in% spike.ins)
##tmp <- tmp[flag,]


ct <- "monocytes"
ct <- "memory.CD4.T.cells"
meth <- "CIBERSORT"
sc <- "fine"
tmp <- subset(res, (measured %in% spike.ins) & (subchallenge == sc))
foo <- subset(tmp, measured > 0)
tmp <- subset(tmp, cell.type %in% unique(foo$cell.type))
flag <- (tmp$measured %in% spike.ins)
tmp <- tmp[flag,]
sub <- subset(tmp, (cell.type == ct) & (subchallenge == sc) & (repo_name == meth) & (dataset.name %in% c("A", "B")) & (measured %in% spike.ins))
## sub$measured <- factor(sub$measured, levels = sort(unique(sub$measured)))
sub <- sub[order(sub$measured),]
sci <- formatC(as.numeric(as.character(sub$measured)), format="e", digits=2)
sub$measured <- factor(sci, levels = unique(sci))
## g <- ggplot(data = sub, aes(x = measured, y = prediction)) + geom_boxplot(aes(group = measured))



## create all plots
## figure out how to subset these results (fine)
## figure out how to subset these results (coarse)



