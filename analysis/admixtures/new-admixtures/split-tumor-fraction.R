library(openxlsx)

## Alternate the 'tumor.fraction' in the admixtures into 'breast' and 'CRC'
prefix <- "060519"
prefix <- "061819"
admixture.file <- paste0(prefix, "-simplified-admixtures.xlsx")
output.admixture.file <- paste0(prefix, "-simplified-admixtures-breast-crc.xlsx")



distribute.tumor.fraction <- function(tbl) {
  tbl$breast <- tbl$tumor.fraction
  tbl$CRC <- 0
  indices <- seq(from=1, to=nrow(tbl), by=2)
  tbl$CRC[indices] <- tbl$breast[indices]
  tbl$breast[indices] <- 0
  tbl <- tbl[, !(colnames(tbl) == "tumor.fraction")]
  cols <- colnames(tbl)
  cols <- c("breast", "CRC", cols[!(cols %in% c("breast", "CRC"))])
  tbl <- tbl[, cols]
  tbl
}


## Sheet 1 has the biological admixtures
biological.admixtures <- read.xlsx(admixture.file, sheet=1)
biological.admixtures <- distribute.tumor.fraction(biological.admixtures)

## Sheet 2 has the random admxitures
random.admixtures <- read.xlsx(admixture.file, sheet=2)
random.admixtures <- distribute.tumor.fraction(random.admixtures)

l <- list("bio" = biological.admixtures, "rand-unconstrained-zero" = random.admixtures)
write.xlsx(l, output.admixture.file)
     