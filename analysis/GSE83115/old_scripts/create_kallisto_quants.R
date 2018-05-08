library(tidyverse)
library(doMC)

#ec2
tmp_dir   <- "/home/ubuntu/tmp/"
kallisto_dir <- "/home/ubuntu/kallisto_cwl/"
fastq_file <- "fastq.tsv"
index_file <- "../GRCH38.idx"    

setwd(tmp_dir)
n_cores = detectCores() - 1




file_df <- read_tsv(fastq_file)
    


run_kallisto_cwl <- function(df){
    cwl_command <- str_c("cwltool", 
                         str_c(kallisto_dir, "quant.cwl"),
                         "--index", index_file,
                         "--threads", n_cores,
                         "--fastq1", df$paths[[1]],
                         "--fastq2", df$paths[[2]], 
                         sep = " ")
    print(cwl_command)
    system(cwl_command)
    cwl_command2 <- str_c("cwltool", 
                          str_c(kallisto_dir, "h5dump.cwl"),
                          "--h5 abundance.h5",
                          sep = " ")
    system(cwl_command2)
    print(cwl_command2)
    abundance_file <- str_c(df$sample_name[[1]], "_kallisto_abundance.tsv")
    file.rename("abundance.tsv", abundance_file)
    walk(c("abundance.h5", "run_info.json"), file.remove)
}

file_df %>% 
    split(.$sample_name) %>% 
    walk(run_kallisto_cwl)











