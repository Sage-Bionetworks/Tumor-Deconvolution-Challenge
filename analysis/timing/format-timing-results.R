suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))

synLogin()

## The predictions are those from the 'rerun-validation' / 'post-competitive' phase, i.e.,
## where the coarse- and fine-grained challenge use the same data.
## This is as opposed to the original competitive phase (against which the
## methods were ranked) where the coarse- and fine-grained challenges
## differed.

predictions.synId <- "syn22314641"

## These are the final predictions from the competitive phase (i.e., when coarse- and fine-grained
## datasets differed)
# predictions.synId <- "syn22314641"


## Read in the "post-competitive" predictions 
obj <- synGet(predictions.synId, downloadFile=TRUE)
res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

source("../utils.R")

res.all <- assign.result.team.names.and.rounds(res.all)

# These are the run times for the cancer-associated Wu and Pelka datasets
# model_run_times_newdataset.csv  

# Instead load the times for the original Challenge runs
timing.tbl <- read.table("model_run_times_original.csv", sep=",", header=TRUE)

flag <- grepl(timing.tbl$team_name, pattern="mcpcounter|timer|quantiseq|xcell|cibersort|epic")
# In res.all, the comparator methods have their objectId set to the team_name
timing.tbl[flag,"id"] <-  timing.tbl[flag,"team_name"]

# coarse is incorrectly and inconsistently spelled as 'course' for some (but not all) methods
# and this is different between timing.tbl and res.all.
# make it consistent (it needn't be correct)
timing.tbl$id <- gsub(timing.tbl$id, pattern="mcpcounter_coarse", replacement="mcpcounter_course")
timing.tbl$id <- gsub(timing.tbl$id, pattern="quantiseq_coarse", replacement="quantiseq_course")
timing.tbl$id <- gsub(timing.tbl$id, pattern="xcell_coarse", replacement="xcell_course")
timing.tbl$id <- gsub(timing.tbl$id, pattern="cibersort_coarse", replacement="cibersort_course")

stopifnot(all(timing.tbl$id %in% res.all$objectId))

res.all <- subset(res.all, submission != "latest")
final.timing.tbl <- merge(timing.tbl[, c("id", "exec_time.s.")], unique(res.all[, c("objectId", "subchallenge", "method.name", "submission")]),
                          by.x = c("id"), by.y = c("objectId"))

stopifnot(nrow(final.timing.tbl) == nrow(timing.tbl))

final.timing.tbl <- final.timing.tbl[, c("method.name", "subchallenge", "submission", "exec_time.s.")]
colnames(final.timing.tbl) <- c("method.name", "subchallenge", "submission", "exec.time.seconds")

# Manually append the run time of CIBERSORTx, which is in 
# more validation-analysis/csx.out | grep Tim
# Time difference of 9.442667 mins
# Use this both for fine and coarse-grained
csx.time <- 9.442667 * 60

final.timing.tbl <- rbind(final.timing.tbl, c("CIBERSORTx", "fine", "1", csx.time), c("CIBERSORTx", "coarse", "1", csx.time))


final.timing.tbl <- final.timing.tbl[order(final.timing.tbl$method.name, final.timing.tbl$subchallenge, final.timing.tbl$submission, decreasing=FALSE),]

write.table(file="method-run-times.tsv", final.timing.tbl, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

final.timing.tbl$rounded <- round(final.timing.tbl$exec.time.seconds)

# These are the results that go in the main text figure
subset(final.timing.tbl, submission=="1" & method.name %in% c("CIBERSORT", "CIBERSORTx", "xCell", "Biogem", "MCP-counter", "mitten_TDC19", "Aginome-XMU", "DA_505", "quanTIseq", "EPIC"))