library(openxlsx)
library(plyr)
library(dplyr)

## Begin Parameters

## prefix <- "060519"
prefix <- "061819"

## Total number of biological admixtures to generate
num.biological.admixtures <- 60

## Total number of random admixtures to generate
num.random.admixtures <- 96 - num.biological.admixtures

## Generate admixtures using each "pool 1" or "pool 2" vials.
## These are assigned in cell-type-to-vial-map
## Pool 1 uses Stem Express, to the extent possible
num.biological.admixtures.using.pool.1 <- floor(num.biological.admixtures/2)
num.biological.admixtures.using.pool.2 <- num.biological.admixtures - num.biological.admixtures.using.pool.1 
num.random.admixtures.using.pool.1 <- floor(num.random.admixtures/2)
num.random.admixtures.using.pool.2 <- num.random.admixtures - num.random.admixtures.using.pool.1 

## Admixture file defines proportions for each cell type.
## Use the cell-type-to-vial map file to assign those proportions to particular vials

admixture.file <- paste0(prefix, "-simplified-admixtures-breast-crc.xlsx")
output.admixture.file <- paste0(prefix, "-simplified-admixtures-breast-crc-vials.xlsx")
cell.type.to.vial.file <- "cell-type-to-vial-map.xlsx"
output.pdf <- paste0(prefix, "-simplified-admixtures-breast-crc-vials.pdf")

## End Parameters

biological.admixture.lengths <-
  list("1" = num.biological.admixtures.using.pool.1,
       "2" = num.biological.admixtures.using.pool.2)

random.admixture.lengths <-
  list("1" = num.random.admixtures.using.pool.1,
       "2" = num.random.admixtures.using.pool.2)

if(sum(unlist(biological.admixture.lengths)) != num.biological.admixtures) {
  stop(paste0("Need to split biological admixtures across pools such that sum ",
              "is total number of desired biological admixtures\n"))
}

if(sum(unlist(random.admixture.lengths)) != num.random.admixtures) {
  stop(paste0("Need to split random admixtures across pools such that sum ",
              "is total number of desired random admixtures\n"))
}

## Sheet 1 has the biological admixtures
biological.admixtures <- read.xlsx(admixture.file, sheet=1)

## Sheet 2 has the random admxitures
random.admixtures <- read.xlsx(admixture.file, sheet=2)

cell.type.to.vial.map <- read.xlsx(cell.type.to.vial.file, sheet=1)
cell.type.to.vial.map <- cell.type.to.vial.map[, c("TASC.vial.num", "Standard.description", "pool", "admixture.col", "vendor")]
cell.type.to.vial.map <- na.omit(cell.type.to.vial.map)
if(length(unique(cell.type.to.vial.map$pool)) != 2) {
  stop(paste0("Was expecting two pools (1,2), but instead have: ",
              paste(unique(cell.type.to.vial.map$pool, collapse=", ")), "\n"))
}
if(!all(unique(cell.type.to.vial.map$pool) == c(1,2))) {
  stop(paste0("Was expecting two pools (1,2), but instead have: ",
              paste(unique(cell.type.to.vial.map$pool, collapse=", ")), "\n"))
}


if(num.biological.admixtures > nrow(biological.admixtures)) {
  stop(paste0("Num requested biological admixtures (", num.biological.admixtures,
              ") exceeds total possible (", nrow(biological.admixtures), ")\n"))
    
}

if(num.random.admixtures > nrow(random.admixtures)) {
  stop(paste0("Num requested random admixtures (", num.random.admixtures,
              ") exceeds total possible (", nrow(random.admixtures), ")\n"))
    
}

assign.admixtures.to.vials.and.pools <- function(tbl, cell.type.to.vial.map, num.admixtures.per.pool) {
  pool.names <- names(num.admixtures.per.pool)
  names(pool.names) <- pool.names
  pool.assignments <- Reduce("c", lapply(pool.names, function(pool) rep(pool, num.admixtures.per.pool[[pool]])))
  tbl <- tbl[1:length(pool.assignments),]
  tbl$pool <- pool.assignments
  rownames(tbl) <- paste0("M", 1:nrow(tbl))
  llply(pool.names,
        .fun = function(pool.name) {
                 sub <- subset(tbl, pool == pool.name)
		 sub <- sub[, !(colnames(sub) == "pool")]
                 sub.map <- subset(cell.type.to.vial.map, pool == pool.name)
		 rownames(sub.map) <- sub.map$admixture.col
		 vials <- sub.map[colnames(sub), "TASC.vial.num"]
		 description <- sub.map[colnames(sub), "Standard.description"]
		 vendor <- sub.map[colnames(sub), "vendor"]		 
		 ret <- rbind(vials, description, vendor, sub)
		 rownames(ret)[1:3] <- c("vial", "description", "vendor")
		 ret
	       })
}

biological.admixtures <-
  assign.admixtures.to.vials.and.pools(biological.admixtures, cell.type.to.vial.map, biological.admixture.lengths)
biological.admixtures <-
  llply(biological.admixtures,
        .fun = function(df) {
	         flag <- !(rownames(df) %in% c("vial", "description", "vendor"))
		 rownames(df)[flag] <- paste0("B", rownames(df)[flag])
		 df
	       })

random.admixtures <-
  assign.admixtures.to.vials.and.pools(random.admixtures, cell.type.to.vial.map, random.admixture.lengths)
random.admixtures <-
  llply(random.admixtures,
        .fun = function(df) {
	         flag <- !(rownames(df) %in% c("vial", "description", "vendor"))
		 rownames(df)[flag] <- paste0("R", rownames(df)[flag])
		 df
	       })

l <- list("bio1" = biological.admixtures[[1]],
          "bio2" = biological.admixtures[[2]],
          "rand1" = random.admixtures[[1]],
          "rand2" = random.admixtures[[2]])
write.xlsx(l, output.admixture.file, row.names=TRUE)

library(ggplot2)
library(reshape2)

nms <- names(l)
names(nms) <- nms
gs <-
  llply(nms,
        .fun = function(nm) {
                 df <- l[[nm]]
		 flag <- !(rownames(df) %in% c("vial", "description", "vendor"))
		 df <- df[flag, ]
		 df <- melt(as.matrix(df))
		 colnames(df) <- c("sample", "cell.type", "fraction")
		 df$fraction <- as.numeric(as.character(df$fraction))
		 g <- ggplot(data = df, aes(x = fraction, cell.type = cell.type))
		 g <- g + geom_histogram()
##                 g <- g + geom_point()
		 g <- g + facet_wrap(~ cell.type)
		 g <- g + ggtitle(nm)
		 dat.text <- ddply(df, .variables = "cell.type", .fun = function(df) data.frame(label = paste0(length(unique(df$fraction)), " unique")))
		 g <- g + geom_text(data = dat.text, mapping = aes(x = 0.25, y = 5, label = label))
		 
		 g
               })

pdf(output.pdf)
l_ply(gs, .fun = function(g) print(g))

flag1 <- !(rownames(biological.admixtures[[2]]) %in% c("vial", "description", "vendor"))
flag2 <- !(rownames(random.admixtures[[2]]) %in% c("vial", "description", "vendor"))
l <- list("bio" = rbind(biological.admixtures[[1]], biological.admixtures[[2]][flag1,]),
          "rand" = rbind(random.admixtures[[1]], random.admixtures[[2]][flag2,]))

nms <- names(l)
names(nms) <- nms
gs <-
  llply(nms,
        .fun = function(nm) {
                 df <- l[[nm]]
		 flag <- !(rownames(df) %in% c("vial", "description", "vendor"))
		 df <- df[flag, ]
		 df <- melt(as.matrix(df))
		 colnames(df) <- c("sample", "cell.type", "fraction")
		 df$fraction <- as.numeric(as.character(df$fraction))
		 g <- ggplot(data = df, aes(x = fraction, cell.type = cell.type))
		 g <- g + geom_histogram()
##                 g <- g + geom_point()
		 g <- g + facet_wrap(~ cell.type)
		 g <- g + ggtitle(nm)
		 dat.text <- ddply(df, .variables = "cell.type", .fun = function(df) data.frame(label = paste0(length(unique(df$fraction)), " unique")))
		 g <- g + geom_text(data = dat.text, mapping = aes(x = 0.25, y = 5, label = label))
		 g
               })

l_ply(gs, .fun = function(g) print(g))
d <- dev.off()
