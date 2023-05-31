suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(ggupset))
#devtools::install_github("krassowski/complex-upset")
suppressPackageStartupMessages(p_load(ComplexUpset))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(cowplot))
suppressPackageStartupMessages(p_load(ComplexHeatmap))
suppressPackageStartupMessages(p_load(xlsx))

set.seed(1234)

source("../utils.R")

synLogin()

method.name.col <- "method.name"
cell.type.col <- "cell.type"
nbins <- 4
na.bin <- nbins + 1

figs.dir <- "figs/"
dir.create(figs.dir, showWarnings = FALSE)

## Read in the healthy and cancer predictions (created by score-cancer-datasets.R)
all.res <- read.table("../cancer-validation/figs/cancer-validation-dataset-all-scores.tsv", sep="\t", header=TRUE)
all.res <- all.res[!(all.res[, method.name.col] %in% c("TIMER", "timer")),]

fill.in.nas <- function(mat, row.name, col.name, value.var) {
  ret <- melt(acast(mat, as.formula(paste0(row.name, "~", col.name)), value.var = value.var))
  colnames(ret) <- c(row.name, col.name, value.var)
  ret
}

# from https://stackoverflow.com/questions/75489458/getting-cutoff-values-for-each-ntile-group
# The best bin value is 1.
bin.value <- function(mat, col.name, higher.is.better = FALSE, nbins = 4) {
  # mat$bin_number <- ntile(factor * mat[, col.name], nbins)
  #na.flag <- is.na(mat[, col.name])
  #mat$bin_number <- NA
  #mat[!na.flag, "bin_number"] <- ntile(rank(factor * mat[!na.flag, col.name], ties.method="min"), nbins)
  qs <- quantile(mat[, col.name], probs = seq(0, 1, length.out = nbins + 1), na.rm=TRUE)
  mat$cut_bin <- cut(mat[, col.name], breaks = qs[!duplicated(qs)], include.lowest = TRUE)
  mat$bin_number <- cut(mat[, col.name], breaks = qs[!duplicated(qs)], include.lowest = TRUE, labels=FALSE)
  if(higher.is.better) {
    mat$bin_number <- max(mat$bin_number, na.rm=TRUE) - mat$bin_number + 1
  }
  return(mat)
}

all.res <-
  ddply(all.res, .variables = c("dataset.name"),
        .fun = function(mat) fill.in.nas(mat, row.name = cell.type.col, col.name = method.name.col, value.var = "cor.p"))

# combine all the challenge validation datasets into a 'Healthy' dataset by taking the mean
flag <- all.res$dataset.name %in% c("Pelka", "Wu")
all.res[!flag,"dataset.name"] <- "Healthy"
#flag <- all.res$dataset.name %in% c("Wu")
#all.res[flag,"dataset.name"] <- "Wu (BRCA)"
#flag <- all.res$dataset.name %in% c("Pelka")
#all.res[flag,"dataset.name"] <- "Pelka (CRC)"
all.res <- 
  ddply(all.res, .variables = c("dataset.name", "method.name", "cell.type"),
        .fun = function(df) {
          data.frame(cor.p = mean(df$cor.p))
        })

# Now take the average and mean over the 3 datasets (healthy, wu, and pelka)
# Note that some datasets do not have 
mean.all.res <- 
  ddply(all.res, .variables = c("method.name", "cell.type"),
        .fun = function(df) {
          data.frame(mean.cor.p = mean(df$cor.p))
        })

mean.all.res <- 
  ddply(mean.all.res,
        .variables = c("cell.type"),
        .fun = function(df) bin.value(df, "mean.cor.p", higher.is.better = TRUE, nbins = nbins))
mean.all.res[is.na(mean.all.res[, "mean.cor.p"]), "bin_number"] <- na.bin

# Cell types that only occur in one dataset won't have an std dev. Remove those
n.datasets.by.cell.type <-
  ddply(unique(all.res[, c("cell.type", "dataset.name")]), .variables = c("cell.type"), 
        .fun = function(df) data.frame(n.datasets = nrow(df)))
cell.types.to.exclude <- as.character(subset(n.datasets.by.cell.type, n.datasets==1)$cell.type)


sd.all.res <- 
  ddply(subset(all.res, !(cell.type %in% cell.types.to.exclude)), .variables = c("method.name", "cell.type"),
        .fun = function(df) {
          data.frame(sd.cor.p = sd(df$cor.p))
        })

sd.all.res <- 
  ddply(sd.all.res,
        .variables = c("cell.type"),
        .fun = function(df) {
            # print(df[1,"cell.type"])
            bin.value(df, "sd.cor.p", higher.is.better = FALSE, nbins = nbins)
        })

sd.all.res[is.na(sd.all.res[, "sd.cor.p"]), "bin_number"] <- na.bin

# Process limit of detection / spike in
# Note: used random for spike-ins in manuscript, separated by coarse- and fine-grained.
# Need to combine. 
# See analysis/in-silico-admixtures/plot-sensitivity-results.R
synIds <- list("sensitivity-analysis-results" = "syn22953003")
tbls <-
  llply(synIds,
        .fun = function(synId) {
          file <- synGet(synId, downloadFile = TRUE)$path
          readRDS(file)
        })
coarse.spikein <- subset(tbls[["sensitivity-analysis-results"]][["1"]]$res.matrices$coarse, mixture.type == "Random")
fine.spikein <- subset(tbls[["sensitivity-analysis-results"]][["1"]]$res.matrices$fine, mixture.type == "Random")
coarse.spikein$lod <- coarse.spikein$min.diff.prop * 100
coarse.spikein$subchallenge <- "coarse"
fine.spikein$lod <- fine.spikein$min.diff.prop * 100
fine.spikein$subchallenge <- "fine"
stopifnot(length(unique(coarse.spikein$dataset.name)) == 1)
stopifnot(length(unique(fine.spikein$dataset.name)) == 1)

coarse.spikein <- fill.in.nas(coarse.spikein, row.name = cell.type.col, col.name = method.name.col, value.var = "lod")
fine.spikein <- fill.in.nas(fine.spikein, row.name = cell.type.col, col.name = method.name.col, value.var = "lod")

all.spikein <- rbind(coarse.spikein, fine.spikein)
all.spikein <- subset(all.spikein, !(method.name %in% c("TIMER", "timer")))
# Merge coarse and fine-grained results by taking the mean
all.spikein <-
  ddply(all.spikein,
        .variables = c("cell.type", "method.name"),
        .fun = function(df) data.frame(lod = mean(df$lod)))

all.spikein <- fill.in.nas(all.spikein, row.name = cell.type.col, col.name = method.name.col, value.var = "lod")

all.spikein <- 
  ddply(all.spikein,
        .variables = c("cell.type"),
        .fun = function(df) {
          print(df[1,"cell.type"])
          bin.value(df, "lod", higher.is.better = FALSE, nbins = nbins)
        })

all.spikein[is.na(all.spikein[, "lod"]), "bin_number"] <- na.bin

# Process  spillover
# specificity-analysis/figs/spillover-summary-round-1.tsv
spillover.res <- read.table("../specificity-analysis/figs/spillover-summary-round-1.tsv", sep="\t", header=TRUE)
spillover.res <- fill.in.nas(spillover.res, row.name = cell.type.col, col.name = method.name.col, value.var = "spillover")
spillover.res <- subset(spillover.res, !(method.name %in% c("TIMER", "timer")))
spillover.res <- 
  ddply(spillover.res,
        .variables = c("cell.type"),
        .fun = function(df) {
          # print(df[1,"cell.type"])
          bin.value(df, "spillover", higher.is.better = FALSE, nbins = nbins)
        })

spillover.res[is.na(spillover.res[, "spillover"]), "bin_number"] <- na.bin

# Process sample-level analyses
sample.level.coarse.res <- read.table("../sample-level-analysis/figs/sample-level-round-1-metric.sum.res-coarse.tsv", sep="\t", header=TRUE)
sample.level.fine.res <- read.table("../sample-level-analysis/figs/sample-level-round-1-metric.sum.res-fine.tsv", sep="\t", header=TRUE)

method.anno <- get.method.annotations()

subchallenge.col <- "subchallenge"
round.col <- "submission"
method.anno.round <- get.round.specific.annotations(method.anno, round = "1")
# Exclude score-based output
method.anno.round <- subset(method.anno.round, Output != "Score")
sample.level.coarse.res <- 
  subset(sample.level.coarse.res, method.name %in% 
           subset(method.anno.round, is.na(subchallenge) | (subchallenge == "coarse"))$method.name)
sample.level.fine.res <- 
  subset(sample.level.fine.res, method.name %in% 
           subset(method.anno.round, is.na(subchallenge) | (subchallenge == "coarse"))$method.name)

sample.level.coarse.bin <- bin.value(sample.level.coarse.res, "Pearson", higher.is.better = TRUE, nbins = nbins)
sample.level.fine.bin <- bin.value(sample.level.fine.res, "Pearson", higher.is.better = TRUE, nbins = nbins)


# Process run times
timing.tbl <- read.table("../timing/method-run-times.tsv", sep="\t", header=TRUE)
timing.tbl <- subset(timing.tbl, !(method.name %in% c("TIMER", "timer")))
timing.coarse.raw <- subset(timing.tbl, (submission == 1) & (subchallenge == "coarse"))
timing.fine.raw <- subset(timing.tbl, (submission == 1) & (subchallenge == "fine"))

timing.coarse.bin <- bin.value(timing.coarse.raw, "exec.time.seconds", higher.is.better = FALSE, nbins = nbins)
timing.fine.bin <- bin.value(timing.fine.raw, "exec.time.seconds", higher.is.better = FALSE, nbins = nbins)


# Assemble all metrics
raw.bins <- 
  list("mean" = mean.all.res,
       "sd" = sd.all.res,
       "lod" = all.spikein,
       "spillover" = spillover.res,
       "sample_coarse" = sample.level.coarse.bin,
       "sample_fine" = sample.level.fine.bin,
       "timing_coarse" = timing.coarse.bin,
       "timing_fine" = timing.fine.bin)

for(nm in names(raw.bins)) {
  write.table(file = paste0(figs.dir, "/", "binned-score-", nm, ".tsv"), raw.bins[[nm]], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

nms <- names(raw.bins)
names(nms) <- nms
bin.mats <-
  llply(nms,
        .fun = function(nm) {
          mat <- raw.bins[[nm]]
          if("cell.type" %in% colnames(mat)) {
            ret <- acast(mat, method.name ~ cell.type, value.var = "bin_number")
            colnames(ret) <- paste0(colnames(ret), "_", nm)
          } else {
            ret <- mat[, c("bin_number"), drop=FALSE]
            colnames(ret)[1] <- nm
            rownames(ret) <- mat$method.name
          }
          ret <- data.frame("method.name" = rownames(ret), ret)
          ret
        })

all.bins <- Reduce(function(...) merge(..., all=TRUE, by="method.name"), bin.mats)
rownames(all.bins) <- all.bins$method.name
all.bins[is.na(all.bins)] <- na.bin
all.bins <- all.bins[, !(colnames(all.bins) == "method.name")]

sorted.method.scores <- sort(rowSums(all.bins))



cell.type.cols <- unlist(lapply(colnames(all.bins), function(str) strsplit(str, split="_")[[1]][1]))
coarse.cell.types <- gsub(coarse.cell.types, pattern="\\.", replacement = " ")
coarse.cell.types <- gsub(coarse.cell.types, pattern=" cells", replacement = "")
coarse.cell.types <- coarse.cell.types[!(grepl(x=coarse.cell.types, pattern="sample|timing"))]
cell.type.cols <- gsub(cell.type.cols, pattern="\\.", replacement = " ")

coarse.cols <- cell.type.cols %in% coarse.cell.types
fine.cols <- !(cell.type.cols %in% coarse.cell.types) & !(grepl(cell.type.cols, pattern="sample|timing"))

coarse.cell.type.cols <- cell.type.cols[coarse.cols]
coarse.bins <- all.bins[, coarse.cols]
coarse.bins <- coarse.bins[, order(coarse.cell.type.cols)]
coarse.analysis.cols <- unlist(lapply(colnames(coarse.bins), function(str) strsplit(str, split="_")[[1]][2]))

fine.cell.type.cols <- cell.type.cols[fine.cols]
fine.bins <- all.bins[, fine.cols]
fine.bins <- fine.bins[, order(fine.cell.type.cols)]
fine.analysis.cols <- unlist(lapply(colnames(fine.bins), function(str) strsplit(str, split="_")[[1]][2]))

all.bins <- cbind(coarse.bins, fine.bins, all.bins[, c("sample_coarse", "sample_fine", "timing_coarse", "timing_fine")])
cell.type.cols <- unlist(lapply(colnames(all.bins), function(str) strsplit(str, split="_")[[1]][1]))
cell.type.cols <- gsub(cell.type.cols, pattern="\\.", replacement = " ")
cell.type.cols[grepl(x=cell.type.cols, pattern="sample")] <- "Sample-level"
cell.type.cols[grepl(x=cell.type.cols, pattern="timing")] <- "Runtime"
analysis.cols <- unlist(lapply(colnames(all.bins), function(str) strsplit(str, split="_")[[1]][2]))

analysis.cols[analysis.cols == "mean"] <- "Mean"
analysis.cols[analysis.cols == "sd"] <- "SD"
analysis.cols[analysis.cols == "lod"] <- "LoD"
analysis.cols[analysis.cols == "spillover"] <- "Spillover"
analysis.cols[analysis.cols == "coarse"] <- "Coarse"
analysis.cols[analysis.cols == "fine"] <- "Fine"

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colors = c(cbbPalette[2:5], "#000000")
library(RColorBrewer)
colors <- c(rev(brewer.pal(n=4, name='RdBu')), "#000000")
names(colors) <- as.character(1:5)
cat.colors <- cbbPalette[2:(1+length(unique(analysis.cols)))]
cat.names <- c("Mean", "SD", "LoD", "Spillover", "Coarse", "Fine")
names(cat.colors) <- cat.names
row.labels <- rownames(all.bins[names(sorted.method.scores),])
flag <- row.labels %in% get.comparators.cap()
# Make the comparator methods bold
row.labels[flag] <- paste0("**", row.labels[flag], "**")
htmp <-
  Heatmap(all.bins[names(sorted.method.scores),], cluster_rows = FALSE, cluster_columns = FALSE, col = colors,
          name = "Bin",
          column_split = factor(cell.type.cols, levels = unique(cell.type.cols)), column_title_rot = 45,
          top_annotation = HeatmapAnnotation(Category=factor(analysis.cols, levels = cat.names), 
                                             annotation_legend_param = list(labels = cat.names),
                                             col = list(Category = cat.colors),
                                             annotation_name_gp= gpar(fontsize = 8)),
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8),
          column_names_gp = grid::gpar(fontsize = 8),
          row_names_gp = grid::gpar(fontsize = 8),
          row_labels = gt_render(row.labels))
png(paste0(figs.dir, "/", "binned-score-heatmap.png"))
draw(htmp, merge_legend=TRUE)
d <- dev.off()

pdf(paste0(figs.dir, "/", "binned-score-heatmap.pdf"))
draw(htmp, merge_legend=TRUE)
d <- dev.off()

cols <- colnames(all.bins)
tot.df <- data.frame("total.score" = as.numeric(sorted.method.scores))
rownames(tot.df) <- names(sorted.method.scores)
all.bins.with.tot.scores <- merge(all.bins, tot.df, by = "row.names")
colnames(all.bins.with.tot.scores)[1] <- "method.name"
all.bins.with.tot.scores <- all.bins.with.tot.scores[order(all.bins.with.tot.scores$total.score, decreasing=FALSE),]

write.table(file = paste0(figs.dir, "/", "binned-scores.tsv"), all.bins.with.tot.scores, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat(paste("Exiting successfully\n"))
q(status=0)

