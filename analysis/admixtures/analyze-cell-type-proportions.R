library(MASS)
library(dirichlet)
library(gridExtra)
library(ggplot2)
library(plyr)
library(grid)
library(tidyr)
library(openxlsx)

add.pseudo.count <- TRUE

cytof.10k.file <- "10KImmunomes.CyTOF_\ PBMC.2018-07-09.csv"
flow.pbmc.10k.file <- "10KImmunomes.Flow\ Cytometry_\ PBMC.2018-07-09.csv"
flow.wb.10k.file <- "10KImmunomes.Flow\ Cytometry_\ Whole Blood.2018-07-09.csv"
labs.10k.file <- "10KImmunomes.Lab\ Tests_\ Blood\ Count.2018-07-09.csv"
tcga.file <- "iAtlas-sup1.xlsx"
tcga.absolute.file <- "Archive/TCGA.Kallisto.cibersort.absolute.byDNAmeth.tsv"
tcga.relative.file <- "Archive/TCGA.Kallisto.cibersort.relative.tsv"

cytof.tbl <- read.table(cytof.10k.file, sep=",", header=TRUE, as.is=TRUE)
flow.pbmc.tbl <- read.table(flow.pbmc.10k.file, sep=",", header=TRUE, as.is=TRUE)
flow.wb.tbl <- read.table(flow.wb.10k.file, sep=",", header=TRUE, as.is=TRUE)
labs.tbl <- read.table(labs.10k.file, sep=",", header=TRUE, as.is=TRUE)
tcga.tbl <- read.xlsx(tcga.file)
tcga.abs.tbl <- read.table(tcga.absolute.file, sep="\t", header=TRUE, as.is=TRUE)
tcga.rel.tbl <- read.table(tcga.relative.file, sep="\t", header=TRUE, as.is=TRUE)

## filter out the normals from TCGA samples.
## normal status is encoded in the sample name
## Sample type is encoded in the fourth part of the sample identifier as listed on attached spreadsheet

## e.g. TCGA.04.1514.01A is type "01" which is solid tumor. TCGA.06.0678.11A is type "11" which is solid  normal.
## Also see https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
normal.codes <- c(10, 11, 12, 13, 14)

filter.tcga.normals <- function(tcga.tbl, normal.codes) {
    flag <- unlist(lapply(tcga.tbl$SampleID,
                          function(str) {
                              splits <- unlist(strsplit(str, split="\\."))
                              if(length(splits) != 4) {
                                  stop(paste0("Expected tcga sample to have 4 components: ", str, "\n"))
                              }
                              ## extract the numeric portion
                              num <- as.numeric(gsub(splits[4], pattern="(\\d+).*", replacement="\\1"))
                              ret <- ifelse(num %in% normal.codes, TRUE, FALSE)
                              ret
                          }))
    tcga.tbl[!flag, ]
}

tcga.abs.tbl <- filter.tcga.normals(tcga.abs.tbl, normal.codes)
tcga.rel.tbl <- filter.tcga.normals(tcga.rel.tbl, normal.codes)

tcga.abs.lueko.tbl <- tcga.abs.tbl[, c("SampleID", "CancerType", "TotalLeukocyte")]
tcga.rel.lueko.tbl <- tcga.rel.tbl[, c("SampleID", "CancerType", "TotalLeukocyte")]

calc.CI <- function(vec, alpha) {
    qs <- quantile(vec, c((1-alpha)/2, 1 - ((1-alpha)/2)))
    as.numeric(qs)
}

normalize.mat <- function(mat) {
    flag <- unlist(apply(mat, 1, function(row) any(is.na(row))))
    mat <- mat[!flag, ]

    flag <- unlist(apply(mat, 1, function(row) !(all(row == 0))))
    mat <- mat[flag, ]
    
    ## Adjust proportions to sum to 100
    sums <- unlist(apply(mat, 1, function(row) sum(row)))
    for(i in 1:nrow(mat)) {
        mat[i, ] <- mat[i, ] / sums[i]
    }
    mat
}

## Add mean, std, and 95% CI rows to the matrix
summarize.mat <- function(mat) {

    mat <- normalize.mat(mat)
    
    means <- unname(unlist(apply(mat, 2, function(col) mean(col, na.rm=TRUE))))
    sds <- unname(unlist(apply(mat, 2, function(col) sd(col, na.rm=TRUE))))
    bounds <- apply(mat, 2, function(col) calc.CI(col, 0.95))
    lb <- unname(bounds[1,])
    ub <- unname(bounds[2,])
    mat <- rbind(mat, "mean" = means, "sd" = sds, "lower.bound.95.CI" = lb, "upper.bound.95.CI" = ub)
    mat
}

mats <- list()
cnts <- list()
worksheet.names <- list()
tops <- list()
prefixes <- list()
descriptions <- list()

## Remove duplicate subject ids
cytof.tbl <- cytof.tbl[!duplicated(cytof.tbl$subject_accession),]
labs.tbl <- labs.tbl[!duplicated(labs.tbl$subject_accession),]
flow.pbmc.tbl <- flow.pbmc.tbl[!duplicated(flow.pbmc.tbl$subject_accession),]
flow.wb.tbl <- flow.wb.tbl[!duplicated(flow.wb.tbl$subject_accession),]

set.seed(1234)

cols <- c("Leukocyte.Fraction", "Stromal.Fraction")
rownames(tcga.tbl) <- tcga.tbl$TCGA.Participant.Barcode
df <- tcga.tbl[, c("TCGA.Study", cols)]
df <- na.omit(df)
df <- df[order(df$TCGA.Study),]
tcga.lfs.all <- df
tcga.lfs.mats <- dlply(df, .variables = c("TCGA.Study"),
                       .fun = function(df.study) {
                           mat <- as.matrix(df.study[, cols])
                           rownames(mat) <- rownames(df.study)
                           mat
                       })
tcga.lfs.mats[["all"]] <- as.matrix(df[, cols])
rownames(tcga.lfs.mats[["all"]]) <- rownames(df)
non.solid.cancers <- c("LAML", "LCML", "DLBC", "GBM", "LGG")
flag <- !(df$TCGA.Study %in% non.solid.cancers)
tcga.lfs.mats[["solid"]] <- as.matrix(df[flag, cols])
rownames(tcga.lfs.mats[["solid"]]) <- rownames(df)[flag]

for(nm in names(tcga.lfs.mats)) {
    key <- paste0("tcga-lfs-", nm)
    mats[[key]] <- tcga.lfs.mats[[nm]]
    worksheet.names[[key]] <- key
    tops[[key]] <- paste0("TCGA Leukocyte/Stromal (", nm, ")")
    descriptions[[key]] <- paste0("TCGA leukocyte/stromal fraction data (", nm, " cancer)")
    if(nm == "all") {
        descriptions[[key]] <- paste0("TCGA leukocyte/stromal fraction data (all cancers)")
    } else if(nm == "solid") {
        descriptions[[key]] <- paste0("TCGA leukocyte/stromal data (solid cancers, i.e., excluding ", paste(non.solid.cancers, collapse=", "))
    }

    prefixes[[key]] <- paste0("tcga-lfs-", nm)
}

cols <- c("B.Cells.Memory", "B.Cells.Naive", "Dendritic.Cells.Activated", "Dendritic.Cells.Resting",
          "Eosinophils", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2",
          "Mast.Cells.Activated", "Mast.Cells.Resting", "Monocytes", "Neutrophils",
          "NK.Cells.Activated", "NK.Cells.Resting", "Plasma.Cells", "T.Cells.CD4.Memory.Activated",
          "T.Cells.CD4.Memory.Resting", "T.Cells.CD4.Naive", "T.Cells.CD8", "T.Cells.Follicular.Helper",
          "T.Cells.gamma.delta", "T.Cells.Regulatory.Tregs", "Neutrophils",
          "Eosinophils", "Mast.Cells", "Dendritic.Cells", "Macrophages")

subpopulations <- list("Macrophages" = c("Macrophages.M0", "Macrophages.M1", "Macrophages.M2"),
                       "B.Cells" = c("B.Cells.Memory", "B.Cells.Naive"),
                       "T.Cells.CD4" = c("T.Cells.CD4.Memory.Activated", "T.Cells.CD4.Memory.Resting", "T.Cells.CD4.Naive"),
                       "NK.cells" = c("NK.Cells.Activated", "NK.Cells.Resting"),
                       "Dendritic.Cells" = c("Dendritic.Cells.Activated", "Dendritic.Cells.Resting"),
                       "Mast.Cells" = c("Mast.Cells.Activated", "Mast.Cells.Resting"))
parental <- names(subpopulations)

cols <- sort(unique(cols))
cols <- cols[!(cols %in% c("Lymphocytes", names(subpopulations)))]
df <- tcga.tbl[, c("TCGA.Study", cols)]
df <- na.omit(df)
df <- df[order(df$TCGA.Study),]

for(p in parental) {
    sum <- unlist(apply(df[, subpopulations[[p]]], 1, function(row) sum(as.numeric(row))))
    new.df <- data.frame(sum = unname(sum))
    colnames(new.df) <- p
    df <- cbind(df, new.df)
}

keep.cols <- cols
keep.cols <- c(keep.cols, names(subpopulations))
keep.cols <- keep.cols[!(keep.cols %in% unname(unlist(subpopulations)))]
keep.cols <- unique(sort(keep.cols))
eps <- 10^-5
flag <- abs(rowSums(df[, keep.cols], na.rm=TRUE) - 1) > eps
if(any(flag)) {
    cat("CIBERSORT columns should sum to one\n")
    print(rowSums(df[flag, keep.cols]))
    stop("")
}

tcga.immune.all <- df[, c("TCGA.Study", keep.cols)]
studies <- unique(df$TCGA.Study)
names(studies) <- studies
tcga.immune.mats <- llply(studies, .fun = function(study) {
    flag <- as.character(df$TCGA.Study) == as.character(study)
    ret <- as.matrix(df[flag, keep.cols])
    rownames(ret) <- rownames(df)[flag]
    ret
})
tcga.immune.mats[["all"]] <- as.matrix(df[, keep.cols])
rownames(tcga.immune.mats[["all"]]) <- rownames(df)
flag <- !(df$TCGA.Study %in% c("LAML", "LCML", "DLBC", "GBM", "LGG"))
tcga.immune.mats[["solid"]] <- as.matrix(df[flag, keep.cols])
rownames(tcga.immune.mats[["solid"]]) <- rownames(df)[flag]

head(rowSums(tcga.immune.mats[["all"]]))

for(nm in names(tcga.immune.mats)) {
    key <- paste0("tcga-immune-", nm)
    mats[[key]] <- tcga.immune.mats[[nm]]
    worksheet.names[[key]] <- key
    tops[[key]] <- paste0("TCGA (", nm, ")")
    descriptions[[key]] <- paste0("TCGA immune population data (", nm, " cancer)")
    if(nm == "all") {
        descriptions[[key]] <- paste0("TCGA immune population data (all cancers)")
    } else if(nm == "solid") {
        descriptions[[key]] <- paste0("TCGA immune population data (solid cancers, i.e., excluding ", paste(non.solid.cancers, collapse=", "))
    }
    prefixes[[key]] <- paste0("tcga-immune-", nm)
}

## Repeat the above for absolute and relative TCGA
summarize.tcga <- function(tcga.tbl, cancer.type.col = "CancerType", sums.to.one = TRUE) {

    cols <- c("B.cells.naive", "B.cells.memory", "Plasma.cells", "T.cells.CD8",
              "T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated",
              "T.cells.follicular.helper", "T.cells.regulatory..Tregs.",
              "T.cells.gamma.delta", "NK.cells.resting", "NK.cells.activated",
              "Monocytes", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2",
              "Dendritic.cells.resting", "Dendritic.cells.activated", "Mast.cells.resting",
              "Mast.cells.activated", "Eosinophils", "Neutrophils")

    subpopulations <- list("Macrophages" = c("Macrophages.M0", "Macrophages.M1", "Macrophages.M2"),
                           "B.Cells" = c("B.cells.memory", "B.cells.naive"),
                           "T.Cells.CD4" = c("T.cells.CD4.memory.activated",
                                             "T.cells.CD4.memory.resting", "T.cells.CD4.naive"),
                           "NK.cells" = c("NK.cells.activated", "NK.cells.resting"),
                           "Dendritic.Cells" = c("Dendritic.cells.activated", "Dendritic.cells.resting"),
                           "Mast.Cells" = c("Mast.cells.activated", "Mast.cells.resting"))
    parental <- names(subpopulations)

    cols <- sort(unique(cols))
    cols <- cols[!(cols %in% c("Lymphocytes", names(subpopulations)))]
    df <- tcga.tbl[, c(cancer.type.col, cols)]
    rownames(df) <- tcga.tbl$SampleID
    df <- na.omit(df)
    df <- df[order(df[, cancer.type.col]),]

    for(p in parental) {
        sum <- unlist(apply(df[, subpopulations[[p]]], 1, function(row) sum(as.numeric(row))))
        new.df <- data.frame(sum = unname(sum))
        colnames(new.df) <- p
        df <- cbind(df, new.df)
    }

    keep.cols <- cols
    keep.cols <- c(keep.cols, names(subpopulations))
    keep.cols <- keep.cols[!(keep.cols %in% unname(unlist(subpopulations)))]
    keep.cols <- unique(sort(keep.cols))
    if(sums.to.one) {
        eps <- 10^-5
        flag <- abs(rowSums(df[, keep.cols], na.rm=TRUE) - 1) > eps
        if(any(flag)) {
            cat("CIBERSORT columns should sum to one\n")
            print(rowSums(df[flag, keep.cols]))
            stop("")
        }
    }

    tcga.immune.all <- df[, c(cancer.type.col, keep.cols)]
    studies <- unique(df[, cancer.type.col])
    names(studies) <- studies
    tcga.immune.mats <- llply(studies, .fun = function(study) {
        flag <- as.character(df[, cancer.type.col]) == as.character(study)
        ret <- as.matrix(df[flag, keep.cols])
        rownames(ret) <- rownames(df)[flag]
        ret
    })
    tcga.immune.mats[["all"]] <- as.matrix(df[, keep.cols])
    flag <- !(df[, cancer.type.col] %in% c("LAML", "LCML", "DLBC", "GBM", "LGG"))
    tcga.immune.mats[["solid"]] <- as.matrix(df[flag, keep.cols])

    head(rowSums(tcga.immune.mats[["all"]]))
    return(tcga.immune.mats)
}

tcga.abs.immune.mats <- summarize.tcga(tcga.abs.tbl, sums.to.one = FALSE)
tcga.rel.immune.mats <- summarize.tcga(tcga.rel.tbl, sums.to.one = TRUE)

for(nm in names(tcga.abs.immune.mats)) {
    key <- paste0("tcga-abs-immune-", nm)
    mats[[key]] <- tcga.abs.immune.mats[[nm]]
    worksheet.names[[key]] <- key
    tops[[key]] <- paste0("TCGA (absolute proportions; ", nm, ")")
    descriptions[[key]] <- paste0("TCGA immune population data (absolute proportions; ", nm, " cancer)")
    if(nm == "all") {
        descriptions[[key]] <- paste0("TCGA immune population data (absolute proportions; all cancers)")
    } else if(nm == "solid") {
        descriptions[[key]] <- paste0("TCGA immune population data (absolute proportions; solid cancers, i.e., excluding ", paste(non.solid.cancers, collapse=", "))
    }
    prefixes[[key]] <- paste0("tcga-abs-immune-", nm)
}

for(nm in names(tcga.rel.immune.mats)) {
    key <- paste0("tcga-rel-immune-", nm)
    mats[[key]] <- tcga.rel.immune.mats[[nm]]
    worksheet.names[[key]] <- key
    tops[[key]] <- paste0("TCGA (relative proportions; ", nm, ")")
    descriptions[[key]] <- paste0("TCGA immune population data (relative proportions; ", nm, " cancer)")
    if(nm == "all") {
        descriptions[[key]] <- paste0("TCGA immune population data (relative proportions; all cancers)")
    } else if(nm == "solid") {
        descriptions[[key]] <- paste0("TCGA immune population data (relative proportions; solid cancers, i.e., excluding ", paste(non.solid.cancers, collapse=", "))
    }
    prefixes[[key]] <- paste0("tcga-rel-immune-", nm)
}


cols <- c("LYMPHOCYTE_percent", "MONOCYTE_percent", "NEUTROPHIL_percent")
labs.mat <- as.matrix(labs.tbl[, cols])
## rownames(labs.mat) <- labs.tbl$study_accession

cols <- c("CD4_T_cells", "CD8_T_cells")
flow.wb.mat <- as.matrix(flow.wb.tbl[, cols])
## rownames(flow.wb.mat) <- flow.wb.tbl$study_accession

## rownames(flow.pbmc.tbl) <- flow.pbmc.tbl$study_accession

## Combine Central_Memory and Effector Memory into Memory
flow.pbmc.tbl$Memory_CD4_T_cells <- as.numeric(flow.pbmc.tbl$Central_memory_CD4_T_cells) + as.numeric(flow.pbmc.tbl$Effector_memory_CD4_T_cells)

flow.pbmc.tbl$Memory_CD8_T_cells <- as.numeric(flow.pbmc.tbl$Central_memory_CD8_T_cells) + as.numeric(flow.pbmc.tbl$Effector_memory_CD8_T_cells)

cols <- c("Tregs", "Naive_CD4_T_cells", "Memory_CD4_T_cells",
          "Naive_CD8_T_cells", "Memory_CD8_T_cells", "NK_cells",
          "Naive_B_cells", "Memory_B_cells")

flow.pbmc.tbl$non.na.pops <- unlist(apply(flow.pbmc.tbl[, cols], 1, function(row) paste(names(row[!is.na(row)]), collapse=", ")))


## Combine Central_Memory and Effector Memory into Memory
cytof.tbl$Memory_CD4_T_cells <- as.numeric(cytof.tbl$Central_Memory_CD4_T_cells) + as.numeric(cytof.tbl$Effector_Memory_CD4_T_cells)

cytof.tbl$Memory_CD8_T_cells <- as.numeric(cytof.tbl$Central_Memory_CD8_T_cells) + as.numeric(cytof.tbl$Effector_Memory_CD8_T_cells)

## Check that subpopulations are similar to their parental populations
parental <- list("Monocytes", "B_cells", "CD4_T_cells", "CD8_T_cells")
subpopulations <- list("Monocytes" = c("CD16_neg_monocytes", "CD16_pos_monocytes"),
                       "B_cells" = c("Memory_B_cells", "Naive_B_cells"),
                       "CD4_T_cells" = c("Memory_CD4_T_cells", "Naive_CD4_T_cells"),
                       "CD8_T_cells" = c("Memory_CD8_T_cells", "Naive_CD8_T_cells"))

for(p in parental) {
    sum <- unlist(apply(cytof.tbl[, subpopulations[[p]]], 1, function(row) sum(as.numeric(row))))
    cat(paste0("Checking that ", p, " = sum of ", paste(subpopulations[[p]], collapse=", "), "\n"))
    df <- data.frame(cytof.tbl[, p], sum)
    colnames(df) <- c(p, "sum")
}

cols <- c("Monocytes", "Tregs", "Naive_CD4_T_cells", "Memory_CD4_T_cells",
          "Naive_CD8_T_cells", "Memory_CD8_T_cells", "NK_cells",
          "Naive_B_cells", "Memory_B_cells")
cytof.mat <- as.matrix(cytof.tbl[, cols])
## rownames(cytof.mat) <- cytof.tbl$study_accession

plot.means.and.variance <- function(means, variance, main) {
    df <- data.frame(mean = means, upper = means + sqrt(variance),
                     lower = means - sqrt(variance),
                     cell.type = names(means))
    df$cell.type <- factor(df$cell.type, levels = df$cell.type)
    g <- ggplot(data = df, aes(x = cell.type, y = mean))
    g <- g + geom_bar(stat="identity", color="black")
    g <- g + geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)
    g <- g + xlab("Cell Type")
    g <- g + ylab("Proportion")
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    g <- g + ggtitle(main)
}

plot.means.and.bounds <- function(means, lower.bounds, upper.bounds, main) {
    df <- data.frame(mean = means, upper = upper.bounds,
                     lower = lower.bounds,
                     cell.type = names(means))
    df$cell.type <- factor(df$cell.type, levels = df$cell.type)
    g <- ggplot(data = df, aes(x = cell.type, y = mean))
    g <- g + geom_bar(stat="identity", color="black")
    g <- g + geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)
    g <- g + geom_text(aes(x = cell.type, y = 100, label = round(mean, digits=2)))
    g <- g + xlab("Cell Type")
    g <- g + ylab("Proportion")
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    g <- g + ylim(c(min(lower.bounds-0.5, 0), max(upper.bounds+0.5, 110)))
    g <- g + ggtitle(main)
    list("g" = g, "df" = df)
}

## Mean of the Dirichlet distribution
## from: https://en.wikipedia.org/wiki/Dirichlet_distribution
dir.mean <- function(alpha) alpha / sum(alpha)

## Variance of the Dirichlet distribution
## https://en.wikipedia.org/wiki/Dirichlet_distribution
dir.var <- function(alpha) {
    a0 <- sum(alpha)
    var <- alpha * (a0 - alpha)
    var <- var / ( (a0 * a0) * (a0 + 1) )
    var
}


lm_corr_eqn <- function(df, method = "pearson", display.r2 = FALSE, display.pval = FALSE){
    m <- lm(y ~ x, df);
    ct <- cor.test(df$x, df$y, method = method)
    estimate <- ct$estimate
    if(display.r2 == TRUE) { estimate <- estimate*estimate }
    pval <- ct$p.value
    eq <- NULL
    if((method == "pearson") && (display.r2 == TRUE)) { 
      if(display.pval) { 
        eq <- substitute(italic(r)^2~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(italic(r)^2~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else if((method == "pearson") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(italic(r)~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(italic(r)~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else if((method == "spearman") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(rho~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(rho~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else {
      stop(paste("lm_corr_eqn does not know how to handle method = ", method,  " display.r2 = ", display.r2, "\n"))
    }
    as.character(as.expression(eq));                 
}


plot.correlation <- function(x, y, labels = NULL, colors = NULL, display.r2 = FALSE, method = "pearson", display.pval = FALSE, ...) {
  df <- data.frame(x = x, y = y)
  if(!is.null(labels)) {
    df$labels <- labels
  }
  g <- NULL
  if(is.null(labels)) {
    g <- ggplot(df, aes(x = x, y = y))
  } else {
    g <- ggplot(df, aes(x = x, y = y, label = labels))
  }
  if(!is.null(colors)) {
    g <- g + geom_point(aes(colour = colors))
  } else {
    g <- g + geom_point()
  }
  if(!is.null(labels)) {
    g <- g + geom_text(vjust = "inward", hjust = "inward")
  }
##  g <- g + theme(legend.position="none")
  g <- g + geom_smooth(data = df, aes(x = x, y = y), method='lm')
  x.min <- min(df$x, na.rm=TRUE)
  x.max <- max(df$x, na.rm=TRUE)
  y.min <- min(df$y, na.rm=TRUE)
  y.max <- max(df$y, na.rm=TRUE)

  ylimits <- NULL
  use.ggplot.2.2.1.limit.code <- TRUE
  if(use.ggplot.2.2.1.limit.code) {
    ylimits <- ggplot_build(g)$layout$panel_ranges[[1]]$y.range
  } else {
    ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range
  }

  ##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 1 * (y.max - y.min), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  if(nrow(df) > 2) {
      lbl <- lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval)
      g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = 0.8 * ylimits[2], label = lbl, parse=TRUE, ...)
  }
  g
}

source("fit-proportions.R")

fit.proportions(cytof.mat, top = "10K Immunomes Cytof", prefix="10kimmunomes-cytof")

fit.proportions(labs.mat, top = "10K Immunomes Labs", prefix="10kimmunomes-labs")

fit.proportions(flow.wb.mat, top = "10K Immunomes Flow (Whole Blood)", prefix="10kimmunomes-flow-wb")

nm <- "10k.cytof"
mats[[nm]] <- cytof.mat
worksheet.names[[nm]] <- nm
tops[[nm]] <- "10K Immunomes Cytof"
descriptions[[nm]] <- "Cytof data from the 10K Immunomes Website (Healthy Individuals)"
prefixes[[nm]] <- "10kimmunomes-cytof"

nm <- "10k.coulter"
mats[[nm]] <- labs.mat
worksheet.names[[nm]] <- nm
tops[[nm]] <- "10K Immunomes Coulter Counter"
descriptions[[nm]] <- "Coulter Counter data from the 10K Immunomes Website (Healthy Individuals)"
prefixes[[nm]] <- "10kimmunomes-coulter"

nm <- "10k.flow.wb"
mats[[nm]] <- flow.wb.mat
worksheet.names[[nm]] <- nm
tops[[nm]] <- "10K Immunomes Flow (Whole Blood)"
descriptions[[nm]] <- "Flow data (Whole Blood) from the 10K Immunomes Website (Healthy Individuals)"
prefixes[[nm]] <- "10kimmunomes-flow-wb"

uniq.pop.subset <- unique(flow.pbmc.tbl$non.na.pops)
uniq.pop.subset <- uniq.pop.subset[grepl(uniq.pop.subset, pattern=", ")]
indx <- 1
for(subset in uniq.pop.subset) {
    cols <- unlist(strsplit(subset, split=",[ ]*"))
    mat <- as.matrix(flow.pbmc.tbl[, cols])
    mat <- na.omit(mat)
    dash.str <- paste(cols, collapse="-")
    fit.proportions(mat,
                    top = "10K Immunomes Flow (PBMC)",
                    prefix=paste0("10kimmunomes-flow-pbmc-", dash.str))


    nm <- paste0("10k.flow.pbmc", ".", indx)
    mats[[nm]] <- mat
    worksheet.names[[nm]] <- nm
    tops[[nm]] <- "10K Immunomes Flow (PBMC)"
    descriptions[[nm]] <- paste0("Flow data (PBMC) from the 10K Immunomes Website (Healthy Individuals) with populations: ", paste(cols, collapse = ", "))
    prefixes[[nm]] <- paste0("10kimmunomes-flow-pbmc-", dash.str)
    indx <- indx + 1
}

## Chevrier paper
## Fig 3A identifies the following clusters as associated with the
## following phenotypes:
## CD4: T-2, T-3, T-6, T-10, T-13, T-15, T-17, T-18
## CD8: T-0, T-1, T-4, T-5, T-7, T-9, T-11, T-12, T-14, T-16, T-19
## These are further indicated in the Fig 5
## The _second_ panel is TAM:
## M-0, M-1, M-2, M-3, M-4, M-5, M-6, M-7, M-8, M-9, M-10, M-11, M-12,
## M-13, M-14, M-15, M-16
## Note that T cell clusters and TAMs were on different panels
library(openxlsx)
chevrier.file <- "chevrier-table-s6.xlsx"
chevrier.tbl <- read.xlsx(chevrier.file, sheet=2)

parental <- list("TAM", "CD4_Tcells", "CD8_Tcells")
subpopulations <- list("TAM" = c("M-0", "M-1", "M-2", "M-3", "M-4", "M-5",
                                 "M-6", "M-7", "M-8", "M-9", "M-10", "M-11",
                                 "M-12", "M-13", "M-14", "M-15", "M-16"),
                       "CD4_Tcells" = c("T-2", "T-3", "T-6", "T-10", "T-13",
                                        "T-15", "T-17", "T-18"),
                       "CD8_Tcells" = c("T-0", "T-1", "T-4", "T-5", "T-7",
                                        "T-9", "T-11", "T-12", "T-14",
                                        "T-16", "T-19"))
for(p in parental) {
    chevrier.tbl[, p] <- unlist(apply(chevrier.tbl[, subpopulations[[p]]],
                                      1, function(row) sum(as.numeric(row))))
}

cols <- c("TAM", "CD4_Tcells", "CD8_Tcells", "Bcell", "Plasma", "DC", "pDC",
          "NK", "Endothelial", "MSC/Pericyte", "Epithelial", "Granulocyte")
tmp <- as.matrix(chevrier.tbl[, cols])

cols <- c("Bcell", "Plasma", "DC", "pDC", "Granulocyte")
tmp <- as.matrix(chevrier.tbl[, cols])

library(synapser)
synLogin()

## Puram Head and Neck (GSE103322)
## Get the annotations from the Puram Head and Neck data
synId <- "syn11990378"
puram.tbl <- read.table(synGet(synId)$path, sep="\t", header=TRUE)

## Exclude myocytes
puram.tbl <- puram.tbl[!(puram.tbl$cell_type == "myocyte"),]

## Clean up some Puram ids
puram.tbl$sample <- unlist(lapply(puram.tbl$sample, function(x) gsub(x, pattern="^HNSCC_(.+)$", replacement="HNSCC\\1")))
puram.tbl$sample <- unlist(lapply(puram.tbl$sample, function(x) gsub(x, pattern="^HNSCC(.+)$", replacement="HN\\1")))
puram.tbl$patient <- unlist(lapply(puram.tbl$sample, function(x) gsub(x, pattern="^([^_]+)_.+$", replacement="\\1")))
flag <- !grepl(x=puram.tbl$patient, pattern="combo")
puram.tbl <- puram.tbl[flag, ]
flag <- grepl(x=puram.tbl$cell_source, pattern="Lymph")
puram.tbl$patient[flag] <- paste0(puram.tbl$patient[flag], "-LN")
flag <- grepl(x=puram.tbl$cell_source, pattern="Primary")
puram.tbl$patient[flag] <- paste0(puram.tbl$patient[flag], "-P")

tmp <- ddply(puram.tbl, .variables = c("patient", "cell_type"),
             .fun = function(df) data.frame(cnt = nrow(df)))
if(FALSE) {
tmp <- ddply(tmp, .variables = c("patient"),
             .fun = function(df) {
                 if(add.pseudo.count) { df$cnt[df$cnt == 0] <- 1 }
                 tot.cnt <- sum(as.numeric(df$cnt))
                 data.frame(cell_type = df$cell_type, tot.cnt = tot.cnt,
                            freq = as.numeric(df$cnt) / tot.cnt)
             })
tmp <- spread_(tmp[, c("patient", "cell_type", "freq")], key="cell_type", value="freq")
tmp[is.na(tmp)] <- 0
}
if(add.pseudo.count) {
    tmp <- spread_(tmp[, c("patient", "cell_type", "cnt")], key="cell_type", value="cnt", fill = 1)
} else {
    tmp <- spread_(tmp[, c("patient", "cell_type", "cnt")], key="cell_type", value="cnt", fill = 0)
}
rownames(tmp) <- tmp$patient
tmp <- tmp[, !(colnames(tmp) == "patient")]
puram.cnts <- as.matrix(tmp)
tmp <- tmp / rowSums(tmp)

mat <- as.matrix(tmp)
non.immune.cols <- c("tumor", "Endothelial", "Fibroblast")
immune.cols <- colnames(mat)[!(colnames(mat) %in% non.immune.cols)]
cat(paste0("Puram immune cols: ", paste(immune.cols, collapse=", "), "\n"))
cat(paste0("Puram non-immune cols: ", paste(non.immune.cols, collapse=", "), "\n"))
mat <- cbind(mat, "Immune.Cells" = unname(unlist(apply(mat[, immune.cols], 1, function(row) sum(as.numeric(row), na.rm=TRUE)))))
puram.cnts <- cbind(puram.cnts, "Immune.Cells" = unname(unlist(apply(puram.cnts[, immune.cols], 1, function(row) sum(as.numeric(row), na.rm=TRUE)))))
flag <- grepl(rownames(mat), pattern="-LN")
mat.ln <- mat[flag, ]
mat.primary <- mat[!flag, ]
cnts.ln.tmp <- puram.cnts[flag,]
cnts.primary.tmp <- puram.cnts[!flag,]
cnts.ln <- rowSums(cnts.ln.tmp)
names(cnts.ln) <- rownames(cnts.ln.tmp)
cnts.primary <- rowSums(cnts.primary.tmp)
names(cnts.primary) <- rownames(cnts.primary.tmp)
fit.proportions(mat.ln, top = "Puram Head and Neck (LN)", prefix="puram-hn-ln", cnts=cnts.ln)
fit.proportions(mat.primary, top = "Puram Head and Neck (Primary)", prefix="puram-hn-primary", cnts=cnts.primary)

nm <- "puram-ln-non-imm"
cnts[[nm]] <- rowSums(cnts.ln.tmp[, c("Immune.Cells", non.immune.cols)])
names(cnts[[nm]]) <- rownames(cnts.ln.tmp[, c("Immune.Cells", non.immune.cols)])
mats[[nm]] <- mat.ln[, c("Immune.Cells", non.immune.cols)]
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Puram Head and Neck (LN; non-immune populations)"
descriptions[[nm]] <- "Puram scRNA-seq Head and Neck (Lymph Nodes): non-immune populations"
prefixes[[nm]] <- "puram-hn-ln-non-imm"

nm <- "puram-ln-imm-abs"
cnts[[nm]] <- rowSums(cnts.ln.tmp[, immune.cols])
names(cnts[[nm]]) <- rownames(cnts.ln.tmp[, immune.cols])
mats[[nm]] <- mat.ln[, immune.cols]
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Puram Head and Neck (LN; absolute, immune populations)"
descriptions[[nm]] <- "Puram scRNA-seq Head and Neck (Lymph Nodes): absolute, immune populations"
prefixes[[nm]] <- "puram-hn-ln-imm-abs"

rel.nm <- "puram-ln-imm-rel"
cnts[[rel.nm]] <- cnts[[nm]]
mats[[rel.nm]] <- mats[[nm]]
worksheet.names[[rel.nm]] <- rel.nm
tops[[rel.nm]] <- tops[[nm]]
tops[[rel.nm]] <- gsub(x = tops[[rel.nm]], pattern = "absolute", replacement = "relative")
descriptions[[rel.nm]] <- descriptions[[nm]]
descriptions[[rel.nm]] <- gsub(x = descriptions[[rel.nm]], pattern = "absolute", replacement = "relative")
prefixes[[rel.nm]] <- prefixes[[nm]]
prefixes[[rel.nm]] <- gsub(x = prefixes[[rel.nm]], pattern = "abs", replacement = "rel")

nm <- "puram-p-non-imm"
cnts[[nm]] <- rowSums(cnts.primary.tmp[, c("Immune.Cells", non.immune.cols)])
names(cnts[[nm]]) <- rownames(cnts.primary.tmp[, c("Immune.Cells", non.immune.cols)])
mats[[nm]] <- mat.primary[, c("Immune.Cells", non.immune.cols)]
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Puram Head and Neck (primary; non-immune populations)"
descriptions[[nm]] <- "Puram scRNA-seq Head and Neck (Primary): non-immune populations"
prefixes[[nm]] <- "puram-hn-primary-non-imm"

nm <- "puram-p-imm-abs"
cnts[[nm]] <- rowSums(cnts.primary.tmp[, immune.cols])
names(cnts[[nm]]) <- rownames(cnts.primary.tmp[, immune.cols])
mats[[nm]] <- mat.primary[, immune.cols]
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Puram Head and Neck (primary; absolute, immune populations)"
descriptions[[nm]] <- "Puram scRNA-seq Head and Neck (Primary): absolute, immune populations"
prefixes[[nm]] <- "puram-hn-primary-imm-abs"

rel.nm <- "puram-p-imm-rel"
cnts[[rel.nm]] <- cnts[[nm]]
mats[[rel.nm]] <- mats[[nm]]
worksheet.names[[rel.nm]] <- rel.nm
tops[[rel.nm]] <- tops[[nm]]
tops[[rel.nm]] <- gsub(x = tops[[rel.nm]], pattern = "absolute", replacement = "relative")
descriptions[[rel.nm]] <- descriptions[[nm]]
descriptions[[rel.nm]] <- gsub(x = descriptions[[rel.nm]], pattern = "absolute", replacement = "relative")
prefixes[[rel.nm]] <- prefixes[[nm]]
prefixes[[rel.nm]] <- gsub(x = prefixes[[rel.nm]], pattern = "abs", replacement = "rel")

## nm <- "puram-p"
## mats[[nm]] <- mat.primary
## worksheet.names[[nm]] <- nm
## tops[[nm]] <- "Puram Head and Neck (Primary)"
## descriptions[[nm]] <- "Puram scRNA-seq Head and Neck (Primary)"
## prefixes[[nm]] <- "puram-hn-primary"


## Tirosh Melanoma (GSE72056)
synId <- "syn12119636"
tirosh.tbl <- read.table(synGet(synId)$path, sep="\t", header=TRUE)

tmp <- ddply(tirosh.tbl, .variables = c("patient", "cell_type"),
             .fun = function(df) data.frame(cnt = nrow(df)))
if(add.pseudo.count) {
    tmp <- spread_(tmp[, c("patient", "cell_type", "cnt")], key="cell_type", value="cnt", fill = 1)
} else {
    tmp <- spread_(tmp[, c("patient", "cell_type", "cnt")], key="cell_type", value="cnt", fill = 0)
}
rownames(tmp) <- tmp$patient
tmp <- tmp[, !(colnames(tmp) == "patient")]
tirosh.cnts <- as.matrix(tmp)
tmp <- tmp / rowSums(tmp)

mat <- as.matrix(tmp)
tirosh.mat <- mat

nm <- "tirosh"
mats[[nm]] <- mat
cnts[[nm]] <- rowSums(tirosh.cnts)
names(cnts[[nm]]) <- rownames(tirosh.cnts)
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Tirosh Melanoma"
descriptions[[nm]] <- "Tirosh scRNA-seq Melanoma data"
prefixes[[nm]] <- "tirosh-melanoma"

fit.proportions(tirosh.mat, top = "Tirosh Melanoma", prefix="tirosh-melanoma", cnts=cnts[["tirosh"]])

non.immune.cols <- c("CAF", "Endo.", "melanoma")
immune.cols <- colnames(mat)[!(colnames(mat) %in% non.immune.cols)]
cat(paste0("Tirosh immune cols: ", paste(immune.cols, collapse=", "), "\n"))
cat(paste0("Tirosh non-immune cols: ", paste(non.immune.cols, collapse=", "), "\n"))
mat <- cbind(mat, "Immune.Cells" = unname(unlist(apply(mat[, immune.cols], 1, function(row) sum(as.numeric(row), na.rm=TRUE)))))
tirosh.cnts <- cbind(tirosh.cnts, "Immune.Cells" = unname(unlist(apply(tirosh.cnts[, immune.cols], 1, function(row) sum(as.numeric(row), na.rm=TRUE)))))

nm <- "tirosh-non-imm"
cnts[[nm]] <- rowSums(tirosh.cnts[, c("Immune.Cells", non.immune.cols)])
names(cnts[[nm]]) <- rownames(tirosh.cnts[, c("Immune.Cells", non.immune.cols)])
mats[[nm]] <- mat[, c("Immune.Cells", non.immune.cols)]
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Tirosh Melanoma (Non-Immune Populations)"
descriptions[[nm]] <- "Tirosh scRNA-seq Melanoma data (Non-immune populations)"
prefixes[[nm]] <- "tirosh-melanoma-non-imm"

nm <- "tirosh-imm-abs"
cnts[[nm]] <- rowSums(tirosh.cnts[, immune.cols])
names(cnts[[nm]]) <- rownames(tirosh.cnts[, immune.cols])
mats[[nm]] <- mat[, immune.cols]
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Tirosh Melanoma (absolute, immune Populations)"
descriptions[[nm]] <- "Tirosh scRNA-seq Melanoma data (absolute, immune populations)"
prefixes[[nm]] <- "tirosh-melanoma-imm-abs"

rel.nm <- "tirosh-imm-rel"
cnts[[rel.nm]] <- cnts[[nm]]
mats[[rel.nm]] <- mats[[nm]]
worksheet.names[[rel.nm]] <- rel.nm
tops[[rel.nm]] <- tops[[nm]]
tops[[rel.nm]] <- gsub(x = tops[[rel.nm]], pattern = "absolute", replacement = "relative")
descriptions[[rel.nm]] <- descriptions[[nm]]
descriptions[[rel.nm]] <- gsub(x = descriptions[[rel.nm]], pattern = "absolute", replacement = "relative")
prefixes[[rel.nm]] <- prefixes[[nm]]
prefixes[[rel.nm]] <- gsub(x = prefixes[[rel.nm]], pattern = "abs", replacement = "rel")



## Li CRC (GSE81861)
## synId <- "syn11898281"
synId <- "syn11898215"  ## this is the expression file
li.tbl <- read.table(synGet(synId)$path, sep="\t", header=TRUE, comment.char="", quote="\"")
li.tbl <- li.tbl[, c("sample_id", "cell_type")]
colnames(li.tbl) <- c("patient", "cell_type")
li.tbl$patient <- gsub(li.tbl$patient, pattern="_Colorectal_Tumor", replacement="")
li.tbl <- na.omit(li.tbl)

tmp <- ddply(li.tbl, .variables = c("patient", "cell_type"),
             .fun = function(df) data.frame(cnt = nrow(df)))
if(add.pseudo.count) {
    tmp <- spread_(tmp[, c("patient", "cell_type", "cnt")], key="cell_type", value="cnt", fill = 1)
} else {
    tmp <- spread_(tmp[, c("patient", "cell_type", "cnt")], key="cell_type", value="cnt", fill = 0)
}
rownames(tmp) <- tmp$patient
tmp <- tmp[, !(colnames(tmp) == "patient")]
li.cnts <- as.matrix(tmp)
tmp <- tmp / rowSums(tmp)

mat <- as.matrix(tmp)
li.mat <- mat

if(FALSE) {
    fit.proportions(mat.ln, top = "Puram Head and Neck (LN)", prefix="puram-hn-ln")
    fit.proportions(mat.primary, top = "Puram Head and Neck (Primary)", prefix="puram-hn-primary")

    fit.proportions(li.mat, top = "Li CRC", prefix="li-crc")
    fit.proportions(tirosh.mat, top = "Tirosh Melanoma", prefix="tirosh-melanoma")
}


nm <- "li"
cnts[[nm]] <- rowSums(li.cnts)
names(cnts[[nm]]) <- rownames(li.cnts)
mats[[nm]] <- mat
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Li CRC"
descriptions[[nm]] <- "Li scRNA-seq CRC data"
prefixes[[nm]] <- "li-crc"

fit.proportions(li.mat, top = "Li CRC", prefix="li-crc", cnts=cnts[["li"]])

non.immune.cols <- c("Endothelial", "Epithelial", "Fibroblast")
immune.cols <- colnames(mat)[!(colnames(mat) %in% non.immune.cols)]
cat(paste0("Li immune cols: ", paste(immune.cols, collapse=", "), "\n"))
cat(paste0("Li non-immune cols: ", paste(non.immune.cols, collapse=", "), "\n"))
mat <- cbind(mat, "Immune.Cells" = unname(unlist(apply(mat[, immune.cols], 1, function(row) sum(as.numeric(row), na.rm=TRUE)))))
li.cnts <- cbind(li.cnts, "Immune.Cells" = unname(unlist(apply(li.cnts[, immune.cols], 1, function(row) sum(as.numeric(row), na.rm=TRUE)))))

nm <- "li-non-imm"
cnts[[nm]] <- rowSums(li.cnts[, c("Immune.Cells", non.immune.cols)])
names(cnts[[nm]]) <- rownames(li.cnts[, c("Immune.Cells", non.immune.cols)])
mats[[nm]] <- mat[, c("Immune.Cells", non.immune.cols)]
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Li CRC (Non-immune populations)"
descriptions[[nm]] <- "Li scRNA-seq CRC data (Non-immune populations)"
prefixes[[nm]] <- "li-crc-non-imm"

nm <- "li-imm-abs"
cnts[[nm]] <- rowSums(li.cnts[, immune.cols])
names(cnts[[nm]]) <- rownames(li.cnts[, immune.cols])
mats[[nm]] <- mat[, immune.cols]
worksheet.names[[nm]] <- nm
tops[[nm]] <- "Li CRC (absolute, immune populations)"
descriptions[[nm]] <- "Li scRNA-seq CRC data (absolute, immune populations)"
prefixes[[nm]] <- "li-crc-imm-abs"

rel.nm <- "li-imm-rel"
cnts[[rel.nm]] <- cnts[[nm]]
mats[[rel.nm]] <- mats[[nm]]
worksheet.names[[rel.nm]] <- rel.nm
tops[[rel.nm]] <- tops[[nm]]
tops[[rel.nm]] <- gsub(x = tops[[rel.nm]], pattern = "absolute", replacement = "relative")
descriptions[[rel.nm]] <- descriptions[[nm]]
descriptions[[rel.nm]] <- gsub(x = descriptions[[rel.nm]], pattern = "absolute", replacement = "relative")
prefixes[[rel.nm]] <- prefixes[[nm]]
prefixes[[rel.nm]] <- gsub(x = prefixes[[rel.nm]], pattern = "abs", replacement = "rel")

nms <- names(mats)
ordered.nms <- c("10k.cytof", "10k.flow.wb", "10k.coulter", sort(nms[grepl(nms, pattern="10k.flow.pbmc")]),
                 "puram-ln-non-imm", "puram-ln-imm-abs", "puram-ln-imm-rel", "puram-p-non-imm", "puram-p-imm-abs", "puram-p-imm-rel", "tirosh-non-imm", "tirosh-imm-abs", "tirosh-imm-rel", "li-non-imm", "li-imm-abs", "li-imm-rel")
## ordered.nms <- c(ordered.nms, "tcga-lfs-all", "tcga-immune-all", "tcga-lfs-solid", "tcga-immune-solid")
ordered.nms <- c(ordered.nms, "tcga-lfs-all", "tcga-abs-immune-all", "tcga-rel-immune-all", "tcga-lfs-solid", "tcga-abs-immune-solid", "tcga-rel-immune-solid")
## tcga.nms <- nms[grepl(nms, pattern="tcga-lfs")]
## tcga.nms <- tcga.nms[!(tcga.nms %in% c("tcga-lfs-all", "tcga-lfs-solid"))]
tcga.nms <- nms[grepl(nms, pattern="tcga-abs")]
tcga.nms <- tcga.nms[grepl(tcga.nms, pattern="COAD") | grepl(tcga.nms, pattern="READ") | grepl(tcga.nms, pattern="BRCA")] 
tcga.nms <- tcga.nms[!(tcga.nms %in% c("tcga-abs-immune-all", "tcga-abs-immune-solid"))]
ordered.nms <- c(ordered.nms, unlist(lapply(tcga.nms, function(str) c(str, gsub(str, pattern="abs", replacement="rel")))))



## Output all of the tables into an excel spreadsheet
## Create a new workbook
wb <- createWorkbook("Admixtures")

## Write out a table of contents
addWorksheet(wb, "contents")
df <- data.frame("tab" = unname(unlist(worksheet.names[ordered.nms])), "description" = unname(unlist(descriptions[ordered.nms])), stringsAsFactors = FALSE)
df <- rbind(df, c("", ""))
df <- rbind(df, c("summary row name", "summary row description"))
df <- rbind(df, c("mean", "subpopulation mean"))
df <- rbind(df, c("sd", "subpopulation standard deviation"))
df <- rbind(df, c("lower.bound.95.CI", "lower bound of 95% confidence interval of subpopulation"))
df <- rbind(df, c("upper.bound.95.CI", "upper bound of 95% confidence interval of subpopulation"))
writeData(wb, sheet = 1, df, colNames = TRUE, rowNames = FALSE)
    
sheet.indx <- 2
for(nm in ordered.nms) {
    print(worksheet.names[[nm]])
    addWorksheet(wb, worksheet.names[[nm]])
    mat <- summarize.mat(mats[[nm]])
    writeData(wb, sheet = sheet.indx, mat, colNames = TRUE, rowNames = TRUE)
    sheet.indx <- sheet.indx + 1
}

saveWorkbook(wb, "admixtures.xlsx", overwrite = TRUE)

library(ReporteRs)
library(ggcorrplot)

hc_cormat_order <- function(cormat, hc.method = "complete") {
  dd <- stats::as.dist((1 - cormat) / 2)
  hc <- stats::hclust(dd, method = hc.method)
  hc$order
}

## Create a ppt presentation
one.pptx.file <- FALSE
if(one.pptx.file) {
    mydoc <- pptx()
    
    mydoc <- addSlide( mydoc, slide.layout = 'Title Slide' )
    mydoc <- addTitle( mydoc, 'Tumor Deconvolution: Admixture Distributions' )
    mydoc <- addSubtitle( mydoc , 'Aug 22, 2018')
    
    mydoc <- addSlide( mydoc, slide.layout = 'Title and Content' )
    mydoc <- addTitle( mydoc, 'Summary')
    mydoc <- addParagraph( mydoc, value = c("Correlations are more often negative under relative than under absolute normalization.", "Correlations (under both relative and absolute normalization) are dependent on leukocyte fraction.", "Log-normal distribution induces large-magnitude, unconstrained proportions."))
}

nms <- ordered.nms
nms <- nms[!grepl(nms, pattern="lfs")]
nms <- nms[!(nms %in% "10k.coulter")]

pptx.pat.prefixes <- list("10k" = "10k", "puram" = "puram", "tirosh" = "tirosh", "li" = "li")
tcga.prefixes <- nms[grepl(x=nms, pattern="tcga-abs")]
tcga.prefixes <- unique(gsub(x = tcga.prefixes, pattern="tcga-abs-", replacement=""))
for(tcga.prefix in tcga.prefixes) {
    pptx.pat.prefixes[[tcga.prefix]] <- paste0("tcga-", tcga.prefix)
}
print(pptx.pat.prefixes)

pptx.prefixes <- rep("", length(nms))
names(pptx.prefixes) <- nms
for(pat in names(pptx.pat.prefixes)) {
    flag <- grepl(x = nms, pattern=pat)
    pptx.prefixes[flag] <- pptx.pat.prefixes[[pat]]
}

stop("stop")

if(any(pptx.prefixes == "")) { print(pptx.prefixes); stop("Some pptx prefix was NULL\n") }

cur.pptx.file <- ""
last.nm <- ""

corr.of.corrs.plots <- list()

for(nm in nms) {

    if(!one.pptx.file && (cur.pptx.file != pptx.prefixes[nm])) {

        if(cur.pptx.file != "") {

            nm1 <- last.nm
            if(nm1 %in% names(corr.of.corrs.plots)) {
                nm2 <- nm1
                if(grepl(x=nm2, pattern="rel")) { nm2 <- gsub(x=nm2, pattern="rel", replacement="abs")
                } else if(grepl(x=nm2, pattern="abs")) { nm2 <- gsub(x=nm2, pattern="abs", replacement="rel")
                } 
                g1 <- corr.of.corrs.plots[[nm1]]
                g2 <- corr.of.corrs.plots[[nm2]]
                top <- tops[[nm1]]
                top <- gsub(x=top, pattern="relative proportions; ", replacement="")
                top <- gsub(x=top, pattern="absolute proportions; ", replacement="")
                title <- paste0(top, " correlation of correlations (as leukocyte frac varies)")
                mydoc <- addSlide( mydoc, slide.layout = "Two Content" )
                mydoc <- addTitle( mydoc, title )
                mydoc <- addPlot( doc = mydoc, fun = plot, x = g1)
                mydoc <- addPlot( doc = mydoc, fun = plot, x = g2)
            }

            
            writeDoc(mydoc, paste0("082218-tumor-deconvolution-", cur.pptx.file, "-admixtures.pptx"))
        }
        cur.pptx.file <- pptx.prefixes[nm]
        mydoc <- pptx()
    
        mydoc <- addSlide( mydoc, slide.layout = 'Title Slide' )
        mydoc <- addTitle( mydoc, paste0("Tumor Deconvolution: Admixture Distributions (", nm, ")") )
        mydoc <- addSubtitle( mydoc , 'Aug 22, 2018')
    }

    
    mat <- na.omit(as.matrix(mats[[nm]]))
    if(ncol(mat) <= 2) { next }

    print(nm)
    ## Only fit proportions if things are normalized to sum to 1.
    ## This is not the case for the TCGA rel proportions
    glist <- NULL
    if(!grepl(nm, pattern="abs")) {
        glist <- fit.proportions(mat, top = tops[[nm]], prefix = prefixes[[nm]], plot.pdf = FALSE, cnts=cnts[[nm]])
        title <- paste0(tops[[nm]], " distributions")
        mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
        mydoc <- addTitle( mydoc, title )
        mydoc <- addPlot( doc = mydoc, fun = plot, x = glist[["fit"]])
    }
    
    drop.cols <- unlist(apply(mat, 2, function(col) all(col == col[1])))
    mat <- mat[, !drop.cols]

    if(!grepl(nm, pattern="abs")) {
        mat <- normalize.mat(mat)
    }
    corr <- round(cor(mat, method = "pearson"), 1)
    ##    g <- ggcorrplot(corr, hc.order = TRUE, type = "lower", lab = TRUE)
    g <- ggcorrplot(corr, hc.order = TRUE, lab = TRUE)
    title <- paste0(tops[[nm]], " correlation heatmap")
##    g <- g + ggtite(title)
    mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
    mydoc <- addTitle( mydoc, title )
    mydoc <- addPlot( doc = mydoc, fun = plot, x = g)

    ## Show both absolute and relative normalization together
    if(grepl(nm, pattern="tcga-abs") || grepl(nm, pattern="puram-ln-imm-abs") ||
       grepl(nm, pattern="puram-p-imm-abs") || grepl(nm, pattern="li-imm-abs") ||
       grepl(nm, pattern="tirosh-imm-abs")) {
        top <- tops[[nm]]
        top <- gsub(x = top, pattern = "absolute", replacement = "absolute and relative")
        title <- paste0(top, " correlation heatmaps")
        mydoc <- addSlide( mydoc, slide.layout = "Two Content" )
        mydoc <- addTitle( mydoc, title )

        g <- g + ggtitle("Absolute")
        mydoc <- addPlot( doc = mydoc, fun = plot, x = g)
        
        nm.rm <- nm
        nm.rm <- gsub(x = nm.rm, pattern="abs", replacement="rel")
        rel.mat <- na.omit(as.matrix(mats[[nm.rm]]))
        if(ncol(rel.mat) <= 2) { stop("Was not expecting relative mat to have < 2 columns\n") }

        rel.drop.cols <- unlist(apply(rel.mat, 2, function(col) all(col == col[1])))
        rel.mat <- rel.mat[, !rel.drop.cols]

        ## Impose order from absolute order
        ordered.cols <- colnames(mat)[hc_cormat_order(corr)]
        rel.mat <- rel.mat[, ordered.cols[ordered.cols %in% colnames(rel.mat)]]


        rel.mat <- normalize.mat(rel.mat)
        rel.corr <- round(cor(rel.mat, method = "pearson"), 1)
        g.rel <- ggcorrplot(rel.corr, hc.order = FALSE, lab = TRUE)

        g.rel <- g.rel + ggtitle("Relative (Imposed Order)")
        mydoc <- addPlot( doc = mydoc, fun = plot, x = g.rel)
    }
    
    ## Now look at correlations based on top and bottom quartiles of purity for tcga
    if(grepl(nm, pattern="tcga-abs") || grepl(nm, pattern="tcga-rel")) {
        leuko.tbl.subset <- subset(tcga.abs.lueko.tbl, SampleID %in% rownames(mat))
        if(grepl(nm, pattern="tcga-rel")) {
            leuko.tbl.subset <- subset(tcga.rel.lueko.tbl, SampleID %in% rownames(mat))
        }
        probs <- seq(from = 0, to = 1, by = 0.2)
        qs <- unname(quantile(leuko.tbl.subset$TotalLeukocyte, na.rm=TRUE, probs = probs))

        ## Plot the disribution of leukocyte infiltration
        ## This will be the same for relative and absolute--so just plot it once
        if(grepl(nm, pattern="tcga-abs")) {
            g <- ggplot(data = leuko.tbl.subset, aes(x = TotalLeukocyte))
            g <- g + geom_density()
            for(q in qs) {
                g <- g + geom_vline(xintercept = q, linetype = "dashed")
            }
            title <- paste0(tops[[nm]], " leukocyte fraction")
            mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
            mydoc <- addTitle( mydoc, title )
            mydoc <- addPlot( doc = mydoc, fun = plot, x = g)
        }
        
        ## Impose clustering order from first quantile
        ordered.cols <- NULL
        corMat <- NULL
        pnames <- unlist(lapply(1:(length(probs)-1), function(i) paste0(probs[i], "-", probs[i+1])))

        intersect.cols <- colnames(mat)
        for(i in 1:(length(probs)-1)) {
            flag <- (leuko.tbl.subset$TotalLeukocyte >= qs[i]) & (leuko.tbl.subset$TotalLeukocyte < qs[i+1])
            mat.sub <- mat[rownames(mat) %in% leuko.tbl.subset$SampleID[flag], ]
##            drop.cols <- unlist(apply(mat.sub, 2, function(col) all(col == 0)))
##            mat.sub <- mat.sub[, !drop.cols]
            cor.mat <- cor(mat.sub, method = "pearson")
            cor.mat[is.na(cor.mat)] <- 0

            ## Here, no columns are dropped
            if(i == 1) {
                ordered.cols <- colnames(cor.mat)[hc_cormat_order(cor.mat)]
                corMat     <- matrix(NA, choose(length(ordered.cols),2), length(probs)-1)
            }
            ## Drop columns here
            mat.sub <- mat[rownames(mat) %in% leuko.tbl.subset$SampleID[flag], ]
            mat.sub <- mat.sub[, ordered.cols]
            drop.cols <- unlist(apply(mat.sub, 2, function(col) all(col == 0)))
            mat.sub <- mat.sub[, !drop.cols]
            intersect.cols <- intersect(intersect.cols, colnames(mat.sub))
        }
        
        for(i in 1:(length(probs)-1)) {
            flag <- (leuko.tbl.subset$TotalLeukocyte >= qs[i]) & (leuko.tbl.subset$TotalLeukocyte < qs[i+1])
            mat.sub <- mat[rownames(mat) %in% leuko.tbl.subset$SampleID[flag], ]
##            drop.cols <- unlist(apply(mat.sub, 2, function(col) all(col == 0)))
##            mat.sub <- mat.sub[, !drop.cols]
            cor.mat <- cor(mat.sub, method = "pearson")
            cor.mat[is.na(cor.mat)] <- 0

            ## Here, no columns are dropped
            if(i == 1) {
                ordered.cols <- colnames(cor.mat)[hc_cormat_order(cor.mat)]
                corMat     <- matrix(NA, choose(length(intersect.cols),2), length(probs)-1)
            }
            print(dim(corMat))
            ## Drop columns here
            mat.sub <- mat[rownames(mat) %in% leuko.tbl.subset$SampleID[flag], ]
            mat.sub <- mat.sub[, ordered.cols]
            drop.cols <- unlist(apply(mat.sub, 2, function(col) all(col == 0)))
            mat.sub <- mat.sub[, !drop.cols]
            cor.mat <- cor(mat.sub, method = "pearson")
            cor.mat[is.na(cor.mat)] <- 0
            ## NB: as.dist excludes diag and upper
            corMat[,i] <- as.vector(as.dist(cor.mat[intersect.cols, intersect.cols]))
            
            corr <- round(cor.mat, 1)
            ##    g <- ggcorrplot(corr, hc.order = TRUE, type = "lower", lab = TRUE)
            ##            g <- ggcorrplot(corr, hc.order = TRUE, lab = TRUE)
            ##    g <- g + ggtite(title)
            title <- paste0(tops[[nm]], " correlation heatmap (", format(qs[i], digits = 2), " <= leuko frac < ", format(qs[i+1], digits = 2), "; n = ", nrow(mat.sub), ")")
            if(i >= 2) {
                mydoc <- addSlide( mydoc, slide.layout = "Two Content" )
            } else {
                mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
            }
            mydoc <- addTitle( mydoc, title )

            g <- ggcorrplot(corr, hc.order = FALSE, lab = TRUE)
            if(i >= 2) {
                g <- g + ggtitle("Imposed Order")
            }
            mydoc <- addPlot( doc = mydoc, fun = plot, x = g)

            ## Also add a plot where we have reclustered
            if(i >= 2) {
                g <- ggcorrplot(corr, hc.order = TRUE, lab = TRUE)
                g <- g + ggtitle("Clustered Order")
                mydoc <- addPlot( doc = mydoc, fun = plot, x = g)
            }
        }

        ## Look at correlations of correlations
        ## save(corMat, file = "corMat.Rd")
        ## cat("Done saving corMat\n")
        ## q(status=0)

        colnames(corMat) <- pnames
        g <- ggcorrplot(cor(corMat), lab = TRUE)
        prop.type <- "Relative"
        if(grepl(x=nm, pattern="abs")) { prop.type <- "Absolute" }
        g <- g + ggtitle(paste0("Correlations of Cell Type ", prop.type, " Fractions\nof Samples with Indicated Leukocyte Fractions"))
        corr.of.corrs.plots[[nm]] <- g
        title <- paste0(tops[[nm]], " correlation of correlations (as leukocyte frac varies)")
        mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
        mydoc <- addTitle( mydoc, title )
        mydoc <- addPlot( doc = mydoc, fun = plot, x = g)
        
    }
    
    if("cor" %in% names(glist)) {
        title <- paste0(tops[[nm]], " correlation of correlations")
        mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
        mydoc <- addTitle( mydoc, title )
        mydoc <- addPlot( doc = mydoc, fun = plot, x = glist[["cor"]])
    }

    last.nm <- nm
    
}

nm1 <- last.nm
if(nm1 %in% names(corr.of.corrs.plots)) {
    nm2 <- nm1
    if(grepl(x=nm2, pattern="rel")) { nm2 <- gsub(x=nm2, pattern="rel", replacement="abs")
    } else if(grepl(x=nm2, pattern="abs")) { nm2 <- gsub(x=nm2, pattern="abs", replacement="rel")
    } 

    g1 <- corr.of.corrs.plots[[nm1]]
    g2 <- corr.of.corrs.plots[[nm2]]
    top <- tops[[nm1]]
    top <- gsub(x=top, pattern="relative proportions; ", replacement="")
    top <- gsub(x=top, pattern="absolute proportions; ", replacement="")
    title <- paste0(top, " correlation of correlations (as leukocyte frac varies)")
    mydoc <- addSlide( mydoc, slide.layout = "Two Content" )
    mydoc <- addTitle( mydoc, title )
    mydoc <- addPlot( doc = mydoc, fun = plot, x = g1)
    mydoc <- addPlot( doc = mydoc, fun = plot, x = g2)
}


if(one.pptx.file) {
    writeDoc(mydoc, "082218-tumor-deconvolution-admixtures.pptx")
} else {
    writeDoc(mydoc, paste0("082218-tumor-deconvolution-", cur.pptx.file, "-admixtures.pptx"))
}

    

## Do binomial experiment.

## Load in the MCP counter populations
synId <- "syn11918430"
cell.sig.gene.tbl <- read.table(synGet(synId)$path, sep="\t", header=TRUE)
cell.sig.gene.tbl <- subset(cell.sig.gene.tbl, Method == "mcpcounter")

## Get the data that Andrew Lamb downsampled
synId <- "syn13888758"
gene.cnts <- read.table(synGet(synId)$path, sep="\t", header=TRUE)
rownames(gene.cnts) <- gene.cnts$Hugo
gene.cnts <- as.matrix(gene.cnts[, !(colnames(gene.cnts) %in% c("Hugo"))])
tot.cnts <- colSums(as.matrix(gene.cnts[, !(colnames(gene.cnts) %in% c("Hugo"))]))
gene.cnts <- sweep(gene.cnts, 2, FUN = "/", unname(tot.cnts))
colSums(gene.cnts)
gene.cnts <- gene.cnts[rownames(gene.cnts) %in% cell.sig.gene.tbl$Hugo, ]

library(reshape2)
gene.cnts.long <- melt(gene.cnts)
colnames(gene.cnts.long) <- c("Hugo", "Sample", "Frac")

gene.cnts.long <- merge(gene.cnts.long, cell.sig.gene.tbl[, c("Hugo", "cell_type")])

metagene.tbl <- ddply(gene.cnts.long, .variables = c("Sample", "cell_type"),
                      .fun = function(df) {
                          data.frame(median.frac = median(as.numeric(df$Frac)))
                      })

metagene.tbl$patient <- gsub(metagene.tbl$Sample, pattern="^([^_]+).+$", replacement="\\1")
## Sanity check that cell profile (median) fraction is highest in corresponding purified proportion
g <- ggplot(data = metagene.tbl)
g <- g + geom_point(aes(x = cell_type, y = median.frac))
g <- g + facet_wrap( ~ Sample)
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g <- g + ylab("Fraction of Reads (Median of genes in corresponding profile)")
## https://www.synapse.org/#!Synapse:syn12649844
pdf("gse64655-frac-mcp-sigs-vs-sample.pdf")
g
d <- dev.off()

## Read in the ground truth
synId <- "syn12650217"
truth.tbl <- read.table(synGet(synId)$path, sep="\t", header=TRUE)
truth.tbl.long <- melt(truth.tbl)
colnames(truth.tbl.long) <- c("patient", "cell_type", "proportion")
truth.tbl.long$patient <- gsub(truth.tbl.long$patient, pattern="_PBMC_d0", replacement="")
truth.tbl.long$cell_type <- gsub(truth.tbl.long$cell_type, pattern="\\.", replacement="_")
truth.tbl.long$Sample <- paste(truth.tbl.long$patient, truth.tbl.long$cell_type, sep="_")
truth.tbl.long <- subset(truth.tbl.long, cell_type %in% c("B_cells", "T_cells", "NK_cells"))

flag <- (grepl(metagene.tbl$Sample, pattern="B_cells") & (metagene.tbl$cell_type=="B lineage"))
flag <- flag | (grepl(metagene.tbl$Sample, pattern="T_cells") & (metagene.tbl$cell_type=="T cells"))
flag <- flag | (grepl(metagene.tbl$Sample, pattern="NK_cells") & (metagene.tbl$cell_type=="NK cells"))
tmp <- merge(metagene.tbl[flag, ], truth.tbl.long, by=c("Sample", "patient"))
tmp$proportion <- as.numeric(tmp$proportion) / 100

## tot.reads <- c(118448724, exp(rev(seq(from=log(31250), to=log(2000000), by=log(2)))))
tot.reads <- exp(rev(seq(from=log(31250), to=log(2000000), by=log(2))))

d_ply(tmp, .variables = c("patient"),
      .fun = function(df.patient) {
          summaries <- ldply(1:nrow(df.patient),
                               .fun = function(i) {
                                   cell <- df.patient$cell_type.y[i]
                                   prob <- as.numeric(df.patient$median.frac[i])
                                   prob <- prob * as.numeric(df.patient$proportion[i])
                                   tmp2 <- ldply(1:length(tot.reads),
                                                 .fun = function(j) {
                                                     size <- as.numeric(tot.reads[j])
                                                     vec <- rbinom(n = 1000, size = floor(size), prob = prob)
                                                     bounds <- calc.CI(vec, 0.95)
                                                     data.frame(size = tot.reads[j], mean = mean(vec), lower = bounds[1], upper = bounds[2])
                                                 })
                                   tmp2$cell_type <- cell
                                   tmp2
                               })
          g <- ggplot(data = summaries, aes(x = log(size), y = mean))
          g <- g + geom_bar(stat="identity", color="black")
          g <- g + facet_wrap( ~ cell_type)
          g <- g + geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)
          g <- g + xlab("Total Number of Reads")
          g <- g + ylab("Expected Number of Reads for Cell Type")
          g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
          g <- g + geom_hline(yintercept = 0, linetype = "dashed")
          df.scale <- data.frame(breaks = log(unique(summaries$size)), labels = unique(summaries$size))
          g <- g + scale_x_continuous(breaks = df.scale$breaks, labels = round(df.scale$labels))
          lower.bounds <- summaries$lower
          upper.bounds <- summaries$upper
          g <- g + ylim(c(min(lower.bounds-0.5, 0), max(upper.bounds+0.5, 110)))
          pt <- unique(df.patient$patient)
          g <- g + ggtitle(paste0("Binomial Analysis: ", pt))
          pdf(paste0("gse64655-", pt, "-binomial.pdf"), width = 14)
          print(g)
          d <- dev.off()
                                  
      })

pdf("gse64655-ground-truth.pdf")
grid.table(truth.tbl, rows=NULL)
d <- dev.off()

q(status=0)

ct <- "CAF"
for(ct in c("CAF", "Endo.", "melanoma")) {
  lm <- mean(log10(mats[["tirosh-non-imm"]][,ct]))
  ls <- sd(log10(mats[["tirosh-non-imm"]][,ct]))
  print(c(lm,ls))
  print(c(10^(lm-ls), 10^(lm+ls)))
  sam <- 10^rnorm(n=5000,mean=lm, sd=ls)
  print(mean(sam))
  print(sd(sam))
}
