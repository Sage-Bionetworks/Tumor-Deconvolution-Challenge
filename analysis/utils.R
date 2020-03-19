firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

get.in.silico.metadata <- function() {

    synId <- "syn21773017"
    metadata.file <- synGet(synId, downloadFile = TRUE)$path
    df <- read.table(metadata.file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    return(df)
    
    l <-
        list(
            "A" = list("mixture.type" = "Random", "tumor.type" = NA, "subchallenge" = "fine"),
            "B" = list("mixture.type" = "Biological", "tumor.type" = NA, "subchallenge" = "fine"),            
            "G" = list("mixture.type" = "Random", "tumor.type" = NA, "subchallenge" = "coarse"),
            "H" = list("mixture.type" = "Biological", "tumor.type" = NA, "subchallenge" = "coarse"),
            "C" = list("mixture.type" = "Random", "tumor.type" = "CRC", "subchallenge" = "fine"),
            "E" = list("mixture.type" = "Biological", "tumor.type" = "CRC", "subchallenge" = "fine"),            
            "D" = list("mixture.type" = "Random", "tumor.type" = "BRCA", "subchallenge" = "fine"),
            "F" = list("mixture.type" = "Biological", "tumor.type" = "BRCA", "subchallenge" = "fine"),            
            "I" = list("mixture.type" = "Random", "tumor.type" = "CRC", "subchallenge" = "coarse"),
            "K" = list("mixture.type" = "Biological", "tumor.type" = "CRC", "subchallenge" = "coarse"),            
            "J" = list("mixture.type" = "Random", "tumor.type" = "BRCA", "subchallenge" = "coarse"),
            "L" = list("mixture.type" = "Biological", "tumor.type" = "BRCA", "subchallenge" = "coarse")
        )
    df <- ldply(l, .fun = function(entry) as.data.frame(entry, stringsAsFactors = FALSE))
    colnames(df)[1] <- "dataset.name"
    df
}

get.validation.metadata <- function() {
    ## Our original admixture specification includes the vendor for each sample
    synId <- "syn21577258"
    obj <- synGet(synId, downloadFile = TRUE)
    bm1 <- read.xlsx(obj$path, sheet = 1)
    bm1 <- bm1[-c(1:3),]
    bm2 <- read.xlsx(obj$path, sheet = 2)
    bm2 <- bm2[-c(1:3),]
    rm1 <- read.xlsx(obj$path, sheet = 3)
    rm1 <- rm1[-c(1:3),]
    rm2 <- read.xlsx(obj$path, sheet = 4)
    rm2 <- rm2[-c(1:3),]
    
    lst <- list(bm1, bm2, rm1, rm2)

    df <-
        do.call(rbind, lapply(lst,
                              function(entry) {
                                  label = ifelse(entry$CRC > 0, "CRC",
                                          ifelse(entry$breast > 0, "BRCA", ""))                              
                                  data.frame(id = entry[,1], tumor.type = label)
                              }))
    df$batch <- "empty"
    flag <- df$id %in% bm1[,1]
    df$batch[flag] <- "BM1"
    flag <- df$id %in% bm2[,1]
    df$batch[flag] <- "BM2"
    flag <- df$id %in% rm1[,1]
    df$batch[flag] <- "RM1"
    flag <- df$id %in% rm2[,1]
    df$batch[flag] <- "RM2"
    
    df$mixture.type <- "empty"
    flag <- ( df$id %in% bm1[,1] ) | ( df$id %in% bm2[,1] )
    ##    df$mixture.type[flag] <- "BM"
    df$mixture.type[flag] <- "Biological"
    flag <- ( df$id %in% rm1[,1] ) | ( df$id %in% rm2[,1] )
    ## df$mixture.type[flag] <- "RM"
    df$mixture.type[flag] <- "Random"
    
    df$dataset <- "empty"
    flag <- ( df$tumor.type == "BRCA" ) & ( df$mixture.type == "Biological" )
    df$dataset[flag] <- "DS1"
    flag <- ( df$tumor.type == "CRC" ) & ( df$mixture.type == "Biological" )
    df$dataset[flag] <- "DS2"
    flag <- ( df$tumor.type == "BRCA" ) & ( df$mixture.type == "Random" )
    df$dataset[flag] <- "DS3"
    flag <- ( df$tumor.type == "CRC" ) & ( df$mixture.type == "Random" )
    df$dataset[flag] <- "DS4"

    df
}

get.validation.ground.truth <- function() {
    ## Our original admixture specification includes the vendor for each sample
    synId <- "syn21577258"
    obj <- synGet(synId, downloadFile = TRUE)
    bm1 <- read.xlsx(obj$path, sheet = 1)
    bm1 <- bm1[-c(1:3),]
    bm2 <- read.xlsx(obj$path, sheet = 2)
    bm2 <- bm2[-c(1:3),]
    rm1 <- read.xlsx(obj$path, sheet = 3)
    rm1 <- rm1[-c(1:3),]
    rm2 <- read.xlsx(obj$path, sheet = 4)
    rm2 <- rm2[-c(1:3),]
    
    lst <- list(bm1, bm2, rm1, rm2)

    df <-
        do.call(rbind, lapply(lst,
                              function(entry) {
                                  rownames(entry) <- entry[,1]
                                  entry <- entry[,-1]
                                  ret <- melt(as.matrix(entry))
                                  colnames(ret) <- c("sample.id", "cell.type", "measured")
                                  ret$sample.id <- as.character(ret$sample.id)
                                  ret$cell.type <- as.character(ret$cell.type)
                                  ret$measured <- as.numeric(as.character(ret$measured))
                                  ret
                              }))

    md <- get.validation.metadata()
    colnames(md) <- c("sample.id", "tumor.type", "batch", "mixture.type", "dataset.name")
    df <- merge(df, md)
    
    df
}

plot.admixtures <- function(mat) {

    hc <- hclust(dist(mat))
    labels <- hc$labels[hc$order]
    if("CRC" %in% labels) { labels <- c("CRC", labels[labels != "CRC"]) }
    if("breast" %in% labels) { labels <- c("breast", labels[labels != "breast"]) }
    if("Breast" %in% labels) { labels <- c("Breast", labels[labels != "Breast"]) }    
        
    mat <- mat[labels,]
    if(length(which(c("CRC", "breast", "Breast") %in% labels)) == 2) {
        o <- order(mat["CRC",,drop=T], mat["breast",,drop=T])
        mat <- mat[, o]
    } else {
        hc <- hclust(dist(t(mat)))
        mat <- mat[,hc$labels[hc$order]]
    }
    
    row.levels <- rownames(mat)
    col.levels <- colnames(mat)

    df <- melt(mat)
    df$Var1 <- factor(df$Var1, levels = row.levels)
    df$Var2 <- factor(df$Var2, levels = col.levels)
    
    g <- ggplot(data = df, aes(x = Var2, y = Var1, fill = value))
    g <- g + geom_tile()
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   text = element_text(size=15))
    g <- g + scale_fill_gradient2("Proportion", limits = c(0,1), low = "red", high = "blue", mid = "white", na.value = "black")
    g <- g + xlab("Sample") + ylab("Cell Type")
    g <- g + theme(axis.text.x = element_blank())
    g
}

plot.cell.type.correlation.heatmap <- function(df, show.corr.text = FALSE, id.var = "modelId", cor.var = "cor.p",
                                               cor.type.label = "Pearson\nCorrelation", digits = 2, limits = c(-1, 1),
                                               pval.var = NULL) {
    orig.df <- df
    df <- df[, c(id.var, "cell.type", cor.var)]
    df[, id.var] <- as.character(df[, id.var])
    df$cell.type <- as.character(df$cell.type)
    
    cell.type.means <- ddply(df, .variables = c("cell.type"),
                             .fun = function(tmp) {
                                 ret <- data.frame(id = "mean", cor = mean(tmp[, cor.var], na.rm=TRUE))
                                 colnames(ret)[1] <- id.var
                                 colnames(ret)[2] <- cor.var                                 
                                 ret
                             })
    cell.type.means <- cell.type.means[order(cell.type.means[, cor.var]),]
    
    method.means <- ddply(df, .variables = c(id.var),
                          .fun = function(tmp) {
                              ret <- data.frame(cell.type = "mean", cor = mean(tmp[, cor.var], na.rm=TRUE))
                              colnames(ret)[2] <- cor.var
                              ret
                          })
    
    method.means <- method.means[order(method.means[, cor.var]),]


    cell.type.levels <- c(cell.type.means$cell.type[cell.type.means$cell.type != "mean"], "mean")
    id.levels <- c("mean", method.means[method.means[, id.var] != "mean", id.var])

    df <- rbind(df, cell.type.means[, c(id.var, "cell.type", cor.var)])
    df <- rbind(df, method.means[, c(id.var, "cell.type", cor.var)])    

    df$cor.label <- formatC(df[, cor.var], format="f", digits=digits)
    if(!is.null(pval.var)) {
        p_load(gtools)
        df <- merge(df, orig.df[, c(id.var, "cell.type", pval.var)], all.x = TRUE)
        flag <- is.na(df[, pval.var])
        df[flag, pval.var] <- 1
        df$cor.label <- paste0(stars.pval(df[, pval.var]), "\n", df$cor.label)
    }
    df$cell.type <- factor(df$cell.type, levels = cell.type.levels)
    df[, id.var] <- factor(df[, id.var], levels = id.levels)
    g <- ggplot(data = df, aes_string(y = id.var, x = "cell.type", fill = cor.var))
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
##    g <- g + scale_fill_gradient2(paste0(cor.type.label, "\nCorrelation"),
##                                  limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
    g <- g + scale_fill_gradient2(cor.type.label, limits = limits,
                                  low = "red", high = "blue", mid = "white", na.value = "black")
    ## g <- g + theme(text = element_text(size=20))
    g
}

limit.matrix.to.protein.coding <- function(mat, use.symbols = TRUE) {

    use.biomaRt <- TRUE
    exclude.mt <- FALSE
    if(use.biomaRt) {
        suppressPackageStartupMessages(p_load(biomaRt))
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        pc.tbl <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description","chromosome_name"),
                        filters = 'biotype', values = c('protein_coding'), mart = ensembl)
        if(exclude.mt) {
            pc.tbl <- subset(pc.tbl, chrosome_name != "MT")
        }
        keys <- as.character(pc.tbl$external_gene_name)
        if(!use.symbols) { keys <- as.character(pc.tbl$ensembl_gene_id) }
        mat <- mat[rownames(mat) %in% keys,]
        mat <- sweep(mat, 2, colSums(mat), `/`) * 10^6
        return(mat)
    }
    if(exclude.mt) { stop("Have not implemented exclude.mt for annotationHub -- but is easy\n") }
    suppressPackageStartupMessages(p_load(AnnotationHub))
    suppressPackageStartupMessages(p_load(ensembldb))
    ah <- AnnotationHub()
    flag <- (ah$species == "Homo sapiens") & (ah$genome == "GRCh38") & (ah$dataprovider == "Ensembl") & (ah$rdataclass == "EnsDb")
    ah2 <- ah[flag, ]
    ## as.data.frame(mcols(ah2))[1:10,c("title"),drop=FALSE]
    edb <- ah2[["AH73881"]]

    ## keytypes(edb)
    ## columns(edb)
    keys <- keys(edb, "GENENAME")
    columns <- c("GENEID", "ENTREZID", "GENEBIOTYPE")
    tbl <- ensembldb::select(edb, keys, columns, keytype = "GENENAME")
    pc.tbl <- subset(tbl, GENEBIOTYPE == "protein_coding")

    keys <- as.character(pc.tbl$GENENAME)
    if(!use.symbols) { keys <- as.character(pc.tbl$GENEID) }
    mat <- mat[rownames(mat) %in% keys,]
    mat <- sweep(mat, 2, colSums(mat), `/`) * 10^6
    mat
}
