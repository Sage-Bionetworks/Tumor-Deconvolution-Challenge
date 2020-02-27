## This is based on gen-admixtures-061819.R (the code used, I believe, to create the in vitro admixtures)

suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(dplyr))

## Get min/max for tumor, leukocyte, non-leukocyte stromal populations
## from iAtlas data

collapse.tumor.min.max <- function(tbls) {

    tbl <- ldply(tbls, .fun = function(df) df)
    colnames(tbl) <- c("tumor", "population", "min", "max")
    
    pan.cancer.tbl <-
        ddply(tbl,
              .variables = c("population"),
              .fun = function(df) {
                  data.frame(min = min(df$min), max = max(df$max))
              })
    pan.cancer.tbl
}

calc.col.min.max <- function(df, min.prob = 0, max.prob = 1) {
    cols <- colnames(df)
    names(cols) <- cols
    ldply(cols,
          .fun = function(col) {
              data.frame(min = as.numeric(quantile(na.omit(df[, col]), probs = c(min.prob))),
                         max = as.numeric(quantile(na.omit(df[, col]), probs = c(max.prob))))
          })
}

clean.up.min.max.params <- function(df) {
    flag <- df$min < 0
    df$min[flag] <- 0
    flag <- df$max > 1
    df$max[flag] <- 1
    df
}

invert.min.max.params <- function(df) {
    pop.col <- c("population", ".id")
    pop.col <- pop.col[pop.col %in% colnames(df)][1]
    pops <- df[, pop.col]
    df <- df[, !(colnames(df) %in% pop.col)]
    df <- t(df)
    colnames(df) <- pops
    df
}


define.biological.model <- function(populations, cache.dir = ".") {

    use.eosinophils <- "Eosinophils" %in% populations
    use.memory.b.cells <- "Memory_B_cells" %in% populations
    
    ## Load the Lambrechts data
    lambrechts.tbl <-
        read.table(paste0(cache.dir, "/MetaData.txt"), sep="\t", header=TRUE, as.is=TRUE)
    cols <- colnames(lambrechts.tbl)
    cols <- cols[1:(length(cols)-1)]
    lambrechts.tbl <- lambrechts.tbl[, 2:ncol(lambrechts.tbl)]
    colnames(lambrechts.tbl) <- cols
    lambrechts.tbl$annotated.name <-
        gsub(lambrechts.tbl$annotated.name, pattern="[ ]+$", replacement="")
    
    lambrechts.pt.sample.tbl <-
        unique(lambrechts.tbl[, c("Sample", "Patientnumber")])

    ## Careful -- there were adjacent normal samples; exclude them
    lambrechts.tbl <- subset(lambrechts.tbl, CellFromTumor == TRUE)

    pt.sample.col <- "Sample"

    tmp <- ddply(lambrechts.tbl,
                 .variables = c(pt.sample.col, "annotated.name"),
                 .fun = function(df) { data.frame(count = nrow(df)) })

    tmp <-
        ddply(tmp, .variables = c(pt.sample.col),
              .fun = function(df) { data.frame(annotated.name = df$annotated.name,
                                               freq = df$count / sum(df$count)) })
    tmp <- acast(tmp, Sample ~ annotated.name)
    tmp[is.na(tmp)] <- 0

    lambrechts.freq <- as.data.frame(tmp)

    translation <-
        list("tumor.fraction" = c("cancer cells pt 1", "cancer cells pt 2",
                                  "cancer cells pt 3", "cancer cells pt 4",
                                  "cancer cells pt 5"),
             "CD4_T_cells" = "CD4+ T cells",
             "CD8_T_cells" = "CD8+ T cells",       
             "Fibroblasts" = c("COL4A2-expressing fibroblasts",
                               "GABARAP-expressing fibroblasts",
                               "lower quality fibroblasts",
                               "normal  lung fibroblasts",
                               "PLA2G2A-expressing fibroblasts",
                               "COL12A1-expressing fibroblasts",
                               "TFPI2-expressing fibroblasts"),
             "Neutrophils" = c("granulocytes"),
             "Endothelial_cells" = c("lower quality endothelial cell",
                                     "lymphatic EC", "normal endothelial cell",
                                     "tumour endothelial cell"),
             "mono.macro" = c("macrophages"),
             "Dendritic_cells" = c("monocyte-derived dendritic cells"),
             "NK_cells" = c("natural killer cells"),
             "Tregs" = c("regulatory T cells"))

    if(use.memory.b.cells) {
        translation[["B_cells"]] = c("follicular B cells", "MALT B cells", "plasma B cells")
    } else {
        translation[["Naive_B_cells"]] = c("follicular B cells", "MALT B cells", "plasma B cells")        
    }

    
    for(entry in names(translation)) {
        cols <- translation[[entry]]
        if(!all(cols %in% colnames(lambrechts.freq))) {
            print(cols[!(cols %in% colnames(lambrechts.freq))])
        }
        lambrechts.freq[, entry] <-
            unlist(apply(lambrechts.freq[, translation[[entry]],drop=F], 1,
                         function(row) sum(row)))
    }

    lambrechts.freq <- lambrechts.freq[, names(translation)]
    lambrechts.freq <- lambrechts.freq / as.numeric(rowSums(lambrechts.freq))

    ## Use all samples (tumor core, middle, and edge)
    if(FALSE) {
        ## Take the sample from each patient with the maximum tumor fraction
        tmp <- lambrechts.freq
        tmp$sample <- rownames(tmp)
        tmp <- merge(tmp, lambrechts.pt.sample.tbl, by.x = "sample", by.y = "Sample")
        tmp <- tmp[order(tmp$tumor.fraction, decreasing=TRUE),]
        tmp <- tmp[!duplicated(tmp$Patientnumber),]
        tmp <- tmp[, !(colnames(tmp) %in% c("sample", "Patientnumber"))]
        lambrechts.freq <- tmp
    }

    li.crc.cnt.df <- load.li.cell.type.counts(cache.dir)

    ## Can not compare immune (CD45+) and non-immune (CD45-), as these
    ## were isolated separately.
    ## In fact, there was further isolation into CD45-/CD90-/CD31- (deplete fibroblasts and endothelial and enrich for malignant)
    ## and CD45+/CD3+ (to enrich for T cells).
    ## Hence, really not use to quantitate relative abundance.
    puram.hn.cnt.df <- load.puram.cell.type.counts(cache.dir)
    ## Subset to Primary samples
    flag <- grepl(rownames(puram.hn.cnt.df), pattern="-P")
    puram.hn.cnt.df <- puram.hn.cnt.df[flag, ]
    puram.hn.cnt.df[is.na(puram.hn.cnt.df)] <- 0

    iatlas.tbl <- load.iAtlas.supp(cache.dir)
    ## The last few columns "Lymphocytes", "Neutrophils", "Eosinophils", "Mast.Cells", "Dendritic.Cells", and "Macrophages"
    ## are an aggregation of earlier columns--note that "Neutrophils" and "Eosinophils" occur multiple times.
    cols.to.keep <- 1:(which(colnames(iatlas.tbl) == "Lymphocytes")-1)

    iatlas.tbl <- iatlas.tbl[, cols.to.keep]
    iatlas.tbl$tumor.fraction = 1 - iatlas.tbl$Stromal.Fraction
    iatlas.tbl$leukocyte.fraction = iatlas.tbl$Leukocyte.Fraction
    iatlas.tbl$non.leukocyte.stromal.fraction <- 1 - (iatlas.tbl$tumor.fraction + iatlas.tbl$leukocyte.fraction)

    ## cibersort columns
    cs.cols <- c(
        "B.Cells.Memory",                         
        "B.Cells.Naive",                          
        "Dendritic.Cells.Activated",              
        "Dendritic.Cells.Resting",                
        "Eosinophils",                            
        "Macrophages.M0",                         
        "Macrophages.M1",                         
        "Macrophages.M2",                         
        "Mast.Cells.Activated",                   
        "Mast.Cells.Resting",                     
        "Monocytes",                              
        "Neutrophils",                            
        "NK.Cells.Activated",                     
        "NK.Cells.Resting",                       
        "Plasma.Cells",                           
        "T.Cells.CD4.Memory.Activated",           
        "T.Cells.CD4.Memory.Resting",             
        "T.Cells.CD4.Naive",                      
        "T.Cells.CD8",                            
        "T.Cells.Follicular.Helper",              
        "T.Cells.gamma.delta",                    
        "T.Cells.Regulatory.Tregs"
    )
    
    ## Define new aggregated cibersort tables
    iatlas.tbl$Memory_CD4_T_cells <- as.numeric(iatlas.tbl$T.Cells.CD4.Memory.Activated) +
        as.numeric(iatlas.tbl$T.Cells.CD4.Memory.Resting)
    iatlas.tbl$NK_cells <- as.numeric(iatlas.tbl$NK.Cells.Activated) +
        as.numeric(iatlas.tbl$NK.Cells.Resting)
    iatlas.tbl$Dendritic_cells <- as.numeric(iatlas.tbl$Dendritic.Cells.Activated) +
        as.numeric(iatlas.tbl$Dendritic.Cells.Resting)
    iatlas.tbl$Macrophages <- as.numeric(iatlas.tbl$Macrophages.M0) +
        as.numeric(iatlas.tbl$Macrophages.M1) +
        as.numeric(iatlas.tbl$Macrophages.M2)
    if(!use.memory.b.cells) {
        iatlas.tbl$B.Cells <- as.numeric(iatlas.tbl$B.Cells.Naive)
    }
    
    old.cs.leuko.cols <- c("Memory_CD4_T_cells", "T.Cells.CD4.Naive", "NK_cells",
                           "T.Cells.Regulatory.Tregs", "Dendritic_cells", "Monocytes", "Macrophages",
                           "Neutrophils", "Eosinophils", "T.Cells.CD8")
    new.cs.leuko.cols <- c("Memory_CD4_T_cells", "Naive_CD4_T_cells", "NK_cells",
                           "Tregs", "Dendritic_cells", "Monocytes", "Macrophages",
                           "Neutrophils", "Eosinophils", "CD8_T_cells")
    if(use.memory.b.cells) {
        old.cs.leuko.cols <- c("B.Cells.Memory", "B.Cells.Naive", old.cs.leuko.cols)
        new.cs.leuko.cols <- c("Memory_B_cells", "Naive_B_cells", new.cs.leuko.cols)
    } else {
        old.cs.leuko.cols <- c("B.Cells", old.cs.leuko.cols)
        new.cs.leuko.cols <- c("Naive_B_cells", new.cs.leuko.cols)

    }
        

    iatlas.tbls <- dlply(iatlas.tbl, .variables = "TCGA.Study", .fun = function(df) df)


    iatlas.leuko.stromal.purity.tbls <-
        llply(iatlas.tbls,
              .fun = function(df) {
                  ret <- data.frame(tumor.fraction = 1 - df$`Stromal.Fraction`, leukocyte.fraction = df$Leukocyte.Fraction)
                  ret$non.leukocyte.stromal.fraction <- 1 - (ret$tumor.fraction + ret$leukocyte.fraction)
                  ret
              })
    
    ## Can not compare immune (CD45+) and non-immune (CD45-), as these
    ## were isolated separately.
    ## Hence, can not compare relative abundance of immune vs non-immune.
    ## In four tumors (Mel58, 67, 72 and 74), we sequenced primarily the immune infiltrates (CD45+ cells) and there were only zero or one malignant cells by this definition.
    ## This seems also to be true for CY75, which doesn't have any tumor cells.
    tirosh.melanoma.cnt.df <- load.tirosh.cell.type.counts(cache.dir)
    
    tirosh.non.immune.cnt.orig.df <- tirosh.melanoma.cnt.df[, c("CAF", "Endo.", "melanoma")]
    exclude <- c("CY58", "CY67", "CY72", "CY74", "CY75")
    tirosh.non.immune.cnt.orig.df <- tirosh.non.immune.cnt.orig.df[!(rownames(tirosh.non.immune.cnt.orig.df) %in% exclude),]
    tirosh.non.immune.cnt.orig.df[is.na(tirosh.non.immune.cnt.orig.df)] <- 0
    
    old.cols <- c("CAF", "Endo.", "melanoma")
    new.cols <- c("Fibroblasts", "Endothelial_cells", "Cancer")
    ##old.cols <- c("CAF", "Endo.")
    ##new.cols <- c("Fibroblasts", "Endothelial_cells")
    tirosh.non.immune.cnt.df <- tirosh.non.immune.cnt.orig.df[, old.cols]
    colnames(tirosh.non.immune.cnt.df) <- new.cols
    
    azizi.cnt.df <- load.azizi.cell.type.counts(cache.dir)
    old.cols <- c("T.cell", "Neutrophil", "DC", "Monocyte", "Macrophage", "NK.cell")
    if(use.eosinophils) {
        old.cols <- c(old.cols, "Mast.cell")
    }
    old.cols <- c(old.cols, "B.cell")
    ## NB: We're going to treat eosinophils as noise.  Though Azizi doesn't assay them,
    ## let's call Mast cells Eosinophils.
    new.cols <- c("T_cells", "Neutrophils", "Dendritic_cells",
                  "Monocytes", "Macrophages", "NK_cells")
    if(use.eosinophils) {
        new.cols <- c(new.cols, "Eosinophils")
    }
    if(use.memory.b.cells) {
        new.cols <- c(new.cols, "B_cells")
    } else {
        new.cols <- c(new.cols, "Naive_B_cells")
    }
    azizi.cnt.df <- azizi.cnt.df[, old.cols]
    colnames(azizi.cnt.df) <- new.cols
    
    cytof.frac.df <- load.cytof.10k.cell.type.fractions(cache.dir)
    cols <- colnames(cytof.frac.df)
    tmp <- cytof.frac.df[, cols]
    tmp <- tmp / rowSums(tmp)
    tmp <- na.omit(tmp)
    cytof.frac.df <- tmp

    min.prob <- pnorm(-2)
    max.prob <- pnorm(2)
    
    iatlas.leuko.stromal.unif.fits <-
        llply(iatlas.leuko.stromal.purity.tbls,
              .fun = function(df) {
                  cols <- colnames(df)
                  names(cols) <- cols
                  df <- na.omit(df)
                  if(nrow(df) == 0) { return(NULL) }
                  ret <- ldply(cols,
                               .fun = function(col) {
                                   data.frame(min = as.numeric(quantile(df[, col], probs=c(min.prob))),
                                              max = as.numeric(quantile(df[, col], probs=c(max.prob))))
                               })
                  colnames(ret) <- c("population", "min", "max")
                  ret
              })
    
    ##iatlas.leuko.stromal.unif.pan.cancer.orig.fits <-
    ##    collapse.tumor.min.max(iatlas.leuko.stromal.unif.fits)
    
    tmp <- ldply(iatlas.leuko.stromal.purity.tbls)
    tmp <- tmp[, !(colnames(tmp) == ".id")]
    iatlas.leuko.stromal.unif.pan.cancer.orig.fits <- calc.col.min.max(tmp, min.prob = min.prob, max.prob = max.prob)
    

    ## Get min/max for leukocyte populations from CIBERSORT
    cibersort.leuko.unif.fits <-
        llply(iatlas.tbls,
              .fun = function(df) {
                  df <- df[, old.cs.leuko.cols]
                  colnames(df) <- new.cs.leuko.cols
                  cols <- colnames(df)
                  names(cols) <- cols
                  df <- na.omit(df)
                  if(nrow(df) == 0) { return(NULL) }
                  ret <- ldply(cols,
                               .fun = function(col) {
                                   data.frame(min = as.numeric(quantile(na.omit(df[, col]), probs = c(min.prob))),
                                              max = as.numeric(quantile(na.omit(df[, col]), probs = c(max.prob))))
                               })
                  colnames(ret) <- c("population", "min", "max")
                  ret
              })
    
    ##cibersort.leuko.unif.pan.cancer.orig.fits <-
    ##    collapse.tumor.min.max(cibersort.leuko.unif.orig.fits)
    
    tmp <- ldply(iatlas.tbls)
    tmp <- tmp[, !(colnames(tmp) == ".id")]
    tmp <- tmp[, old.cs.leuko.cols]
    colnames(tmp) <- new.cs.leuko.cols
    cibersort.leuko.unif.pan.cancer.orig.fits <- calc.col.min.max(tmp, min.prob = min.prob, max.prob = max.prob)
    
    lambrechts.unif.orig.fits <- calc.col.min.max(lambrechts.freq, min.prob = min.prob, max.prob = max.prob)
    
    tirosh.non.immune.frac.df <- tirosh.non.immune.cnt.df /
        rowSums(tirosh.non.immune.cnt.df)
    rowSums(tirosh.non.immune.frac.df)
    
    tirosh.non.immune.unif.orig.fits <- calc.col.min.max(tirosh.non.immune.frac.df, min.prob = min.prob,
                                                         max.prob = max.prob)
    colnames(tirosh.non.immune.unif.orig.fits) <- c("population", "min", "max")
    
    azizi.leuko.frac.df <-
        azizi.cnt.df / rowSums(azizi.cnt.df)
    
    mm.cols <- c("Monocytes", "Macrophages")
    azizi.mono.macro.frac.df <-
        azizi.cnt.df[, mm.cols] / rowSums(azizi.cnt.df[, mm.cols])
    
    azizi.leuko.unif.orig.fits <- calc.col.min.max(azizi.leuko.frac.df, min.prob = min.prob, max.prob = max.prob)
    colnames(azizi.leuko.unif.orig.fits) <- c("population", "min", "max")
    
    azizi.mono.macro.unif.orig.fits <- calc.col.min.max(azizi.mono.macro.frac.df, min.prob = min.prob, max.prob = max.prob)
    colnames(azizi.mono.macro.unif.orig.fits) <- c("population", "min", "max")
    
    cytof.leuko.unif.orig.fits <- calc.col.min.max(cytof.frac.df, min.prob = min.prob, max.prob = max.prob)
    colnames(cytof.leuko.unif.orig.fits) <- c("population", "min", "max")
    
    cols <- c("Tregs", "Naive_CD4_T_cells", "Memory_CD4_T_cells",
              "Naive_CD8_T_cells", "Memory_CD8_T_cells")
    tmp <- cytof.frac.df[, cols]
    cytof.t.cell.unif.orig.fits <- calc.col.min.max(na.omit(tmp / rowSums(tmp)), min.prob = min.prob,
                                                    max.prob = max.prob)
    colnames(cytof.t.cell.unif.orig.fits) <- c("population", "min", "max")
    
    cols <- c("Naive_CD8_T_cells", "Memory_CD8_T_cells")
    tmp <- cytof.frac.df[, cols]
    cytof.cd8.t.cell.unif.orig.fits <- calc.col.min.max(na.omit(tmp / rowSums(tmp)), min.prob = min.prob,
                                                        max.prob = max.prob)
    colnames(cytof.cd8.t.cell.unif.orig.fits) <- c("population", "min", "max")
    
    cols <- c("Naive_CD4_T_cells", "Memory_CD4_T_cells")
    tmp <- cytof.frac.df[, cols]
    cytof.cd4.t.cell.unif.orig.fits <- calc.col.min.max(na.omit(tmp / rowSums(tmp)), min.prob = min.prob,
                                                        max.prob = max.prob)
    colnames(cytof.cd4.t.cell.unif.orig.fits) <- c("population", "min", "max")
    
    cols <- c("Naive_B_cells", "Memory_B_cells")
    tmp <- cytof.frac.df[, cols]
    cytof.b.cell.unif.orig.fits <- calc.col.min.max(na.omit(tmp / rowSums(tmp)), min.prob = min.prob,
                                                    max.prob = max.prob)
    colnames(cytof.b.cell.unif.orig.fits) <- c("population", "min", "max")
    
    lambrechts.unif.fits <- clean.up.min.max.params(lambrechts.unif.orig.fits)
    azizi.leuko.unif.fits <- clean.up.min.max.params(azizi.leuko.unif.orig.fits)
    azizi.mono.macro.unif.fits <- clean.up.min.max.params(azizi.mono.macro.unif.orig.fits)
    cytof.t.cell.unif.fits <- clean.up.min.max.params(cytof.t.cell.unif.orig.fits)
    cytof.cd8.t.cell.unif.fits <- clean.up.min.max.params(cytof.cd8.t.cell.unif.orig.fits)
    cytof.cd4.t.cell.unif.fits <- clean.up.min.max.params(cytof.cd4.t.cell.unif.orig.fits)
    cytof.b.cell.unif.fits <- clean.up.min.max.params(cytof.b.cell.unif.orig.fits)
    cytof.leuko.unif.fits <- clean.up.min.max.params(cytof.leuko.unif.orig.fits)
    tirosh.non.immune.unif.fits <- clean.up.min.max.params(tirosh.non.immune.unif.orig.fits)
    cibersort.leuko.unif.pan.cancer.fits  <- clean.up.min.max.params(cibersort.leuko.unif.pan.cancer.orig.fits )
    iatlas.leuko.stromal.unif.pan.cancer.fits <- clean.up.min.max.params(iatlas.leuko.stromal.unif.pan.cancer.orig.fits)
    
    lambrechts.unif.fits <- invert.min.max.params(lambrechts.unif.fits)
    azizi.leuko.unif.fits <- invert.min.max.params(azizi.leuko.unif.fits)
    azizi.mono.macro.unif.fits <- invert.min.max.params(azizi.mono.macro.unif.fits)
    cytof.t.cell.unif.fits <- invert.min.max.params(cytof.t.cell.unif.fits)
    cytof.cd8.t.cell.unif.fits <- invert.min.max.params(cytof.cd8.t.cell.unif.fits)
    cytof.cd4.t.cell.unif.fits <- invert.min.max.params(cytof.cd4.t.cell.unif.fits)
    cytof.b.cell.unif.fits <- invert.min.max.params(cytof.b.cell.unif.fits)
    cytof.leuko.unif.fits <- invert.min.max.params(cytof.leuko.unif.fits)
    tirosh.non.immune.unif.fits <- invert.min.max.params(tirosh.non.immune.unif.fits)
    cibersort.leuko.unif.pan.cancer.fits  <- invert.min.max.params(cibersort.leuko.unif.pan.cancer.fits )
    iatlas.leuko.stromal.unif.pan.cancer.fits <- invert.min.max.params(iatlas.leuko.stromal.unif.pan.cancer.fits)
    
    
    ## "Lambrechts" hierarchical model
    ## - Lambrechts -> tumor, Fibroblasts, Endothelial, NK, Tregs, mono.macro, DC, Neutrophils, CD4, CD8, B 
    ## - mono.macro -> Azizi -> mono, macro
    ## - CD4 -> 10k -> CD4 memory / naive
    ## - CD8 -> 10k -> CD8 memory / naive
    ## - B -> 10k -> memory / naive
    
    ## "Azizi" hierarchical model
    ## - iAtlas -> purity, leuko, non-leuko stromal
    ## - non-leuko stromal -> Tirosh -> endo, fibro
    
    ## - leuko -> Azizi -> T, B, neutro, DC, mono, macro, NK, mast
    ## - T -> 10k -> CD4 memory / naive, CD8 memory / naive, Treg
    ## - B -> 10k -> memory / naive
    ## - mast -> eosinophils
    
    azizi.leuko.children <-
        c("T_cells", "B_cells", "Neutrophils", "Dendritic_cells",
          "Monocytes", "Macrophages", "NK_cells")
    if(use.eosinophils) {
        azizi.leuko.children <- c(azizi.leuko.children, "Eosinophils")
    }  
    
    azizi.mono.macro.children <- c("Monocytes", "Macrophages")

    lambrechts.root.children <- 
        c("tumor.fraction", "Fibroblasts", "Endothelial_cells", "Dendritic_cells",
          "NK_cells", "Tregs", "Neutrophils", "CD4_T_cells", "CD8_T_cells", "mono.macro")
    if(use.memory.b.cells) {
        lambrechts.root.children <- c(lambrechts.root.children, "B_cells")
    } else {
        lambrechts.root.children <- c(lambrechts.root.children, "Naive_B_cells")        
    }
    lambrechts.hierarchical.model <-
        list(
            list(parent = "root",
                 children = lambrechts.root.children,
                 dist = "unif", params = lambrechts.unif.fits),
            list(parent = "CD4_T_cells",
                 children = c("Memory_CD4_T_cells", "Naive_CD4_T_cells"),
                 dist = "unif", params = cytof.cd4.t.cell.unif.fits),
            list(parent = "CD8_T_cells",
                 children = c("Memory_CD8_T_cells", "Naive_CD8_T_cells"),
                 dist = "unif", params = cytof.cd8.t.cell.unif.fits),
            list(parent = "mono.macro",
                 children = c("Monocytes", "Macrophages"),
                 dist = "unif", params = azizi.mono.macro.unif.fits))
    if(use.memory.b.cells) {
        lambrechts.hierarchical.model[[length(lambrechts.hierarchical.model)+1]] <-
            list(parent = "B_cells",
                 children = c("Memory_B_cells", "Naive_B_cells"),
                 dist = "unif", params = cytof.b.cell.unif.fits)
    }
    
    azizi.hierarchical.model <-
        list(
            list(parent = "root",
                 children = c("tumor.fraction", "leukocyte.fraction", "non.leukocyte.stromal.fraction"),
                 dist = "unif", params = iatlas.leuko.stromal.unif.pan.cancer.fits),
            list(parent = "non.leukocyte.stromal.fraction",
                 children = c("Fibroblasts", "Endothelial_cells"),
                 dist = "unif", params = tirosh.non.immune.unif.fits),
            list(parent = "leukocyte.fraction",
                 children = azizi.leuko.children,
                 dist = "unif", params = azizi.leuko.unif.fits),
            list(parent = "T_cells",
                 children = c("Memory_CD4_T_cells", "Naive_CD4_T_cells",
                              "Memory_CD8_T_cells", "Naive_CD8_T_cells", "Tregs"),
                 dist = "unif", params = cytof.t.cell.unif.fits)
        )
    if(use.memory.b.cells) {
        azizi.hierarchical.model[[length(azizi.hierarchical.model) + 1]] <-
            list(parent = "B_cells",
                 children = c("Memory_B_cells", "Naive_B_cells"),
                 dist = "unif", params = cytof.b.cell.unif.fits)
    }
    
    ## "Cibersort" hierarchical model
    ## - iAtlas -> purity, leuko, non-leuko stromal
    ## - non-leuko stromal -> Tirosh -> endo, fibro
    
    ## - leuko -> cibersort -> CD4 memory / naive, CD8, Treg, B memory / naive, neutro, DC, mono, macro, NK, eos
    ## - CD8 -> 10k -> CD8 memory / naive

    ## working
    
    cibersort.leuko.children <-
        c("Memory_CD4_T_cells",
          "Naive_CD4_T_cells", "NK_cells", "Tregs", "Dendritic_cells",
          "Monocytes", "Macrophages", "Neutrophils", "CD8_T_cells")
    if(use.memory.b.cells) {
        cibersort.leuko.children <- c(cibersort.leuko.children, "Memory_B_cells", "Naive_B_cells")
    } else {
        cibersort.leuko.children <- c(cibersort.leuko.children, "Naive_B_cells")
    }
    if(use.eosinophils) {
        cibersort.leuko.children <- c(cibersort.leuko.children, "Eosinophils")
    }
    
    cibersort.hierarchical.model <-
        list(
            list(parent = "root",
                 children = c("tumor.fraction", "leukocyte.fraction", "non.leukocyte.stromal.fraction"),
                 dist = "unif", params = iatlas.leuko.stromal.unif.pan.cancer.fits),
            list(parent = "non.leukocyte.stromal.fraction",
                 children = c("Fibroblasts", "Endothelial_cells"),
                 dist = "unif", params = tirosh.non.immune.unif.fits),
            list(parent = "leukocyte.fraction",
                 children = cibersort.leuko.children,
                 dist = "unif", params = cibersort.leuko.unif.pan.cancer.fits),
            list(parent = "CD8_T_cells",
                 children = c("Memory_CD8_T_cells", "Naive_CD8_T_cells"),
                 dist = "unif", params = cytof.cd8.t.cell.unif.fits)
        )

    final.pops <- c(populations, "tumor.fraction")

    azizi.flat.model <- flatten.hierarchical.unif.model(azizi.hierarchical.model, pops = final.pops)
    cibersort.flat.model <- flatten.hierarchical.unif.model(cibersort.hierarchical.model, pops = final.pops)
    lambrechts.flat.model <- flatten.hierarchical.unif.model(lambrechts.hierarchical.model, pops = final.pops)
    
    tmp <- melt(as.matrix(azizi.flat.model))
    tmp <- rbind(tmp, melt(as.matrix(cibersort.flat.model)))
    colnames(tmp) <- c("row", "col", "value")
    
    tmp <- ddply(tmp, .variables = c("row", "col"),
                 .fun = function(df) {
                     if(df$row[1] == "min") {
                         return(df[which.min(df$value),,drop=F])
                     }
                     df[which.max(df$value),,drop=F]                 
                 })
    
    flat.model <- acast(tmp, row ~ col)
    flat.model["min",] <- unlist(lapply(flat.model["min",], function(x) max(x, 0.01)))
    flat.model["max",] <- unlist(lapply(flat.model["max",], function(x) max(x, 0.01)))
    non.tumor.cols <- colnames(flat.model)[!(colnames(flat.model) == "tumor.fraction")]
    flat.model["max",non.tumor.cols] <- unlist(lapply(flat.model["max",non.tumor.cols], function(x) min(x, 0.5)))
    flat.model
}




