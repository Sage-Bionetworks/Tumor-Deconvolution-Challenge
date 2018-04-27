require(tidyverse)
require(magrittr)

# Requires fastq-tools to be instaled and added to the env path to be called by
# this script:
# https://github.com/dcjones/fastq-tools

# All files will be written to the current working directory

## df
# must be a data_frame with the columns:
# -fraction
# -p1_fastq_file
# -p2_fastq_file 

# df <- data_frame(
#     "fraction" = c(.2, .3, .5),
#     "p1_fastq_file" = 
#         c("sample1_pair1.fastq", 
#           "sample2_pair1.fastq",
#           "sample3_pair1.fastq"),
#     "p2_fastq_file" = 
#         c("sample1_pair2.fastq", 
#           "sample2_pair2.fastq",
#           "sample3_pair2.fastq"))

# each row of the dataframe must come from one sample
# fraction: fraction of reads in the final fastq that will come from this sample
# p1_fastq_file: path to the first fastq file for this sample
# p2_fastq_file: path to the second fastq file for this sample

## other_args 
# must be a string of other args that will get sent to fastq-sample
# https://homes.cs.washington.edu/~dcjones/fastq-tools/fastq-sample.html

combine_paired_fastq_files <- function(df, other_args = ""){
    
    df2 <- df %>% 
        mutate(n_reads = map_int(p1_fastq_file, find_fastq_n_reads)) %>% 
        mutate(n_reads2 = map_int(p2_fastq_file, find_fastq_n_reads)) %>% 
        mutate(n_samples = as.integer(mean(n_reads) * fraction)) %>% 
        mutate(prefix = str_c("tmp", 1:n())) %>% 
        mutate(p1_sample_file = str_c(prefix, ".1.fastq")) %>% 
        mutate(p2_sample_file = str_c(prefix, ".2.fastq")) %>% 
        mutate(sample_command = create_fastq_sample_commands(., other_args))
    if(any(df2$n_reads != df2$nreads2)) stop("reads per pair don't match")
    walk(df2$sample_command, system)
    create_fastq_merge_command(df2$p1_sample_file, "p1") %>%
        system
    create_fastq_merge_command(df2$p1_sample_file, "p2") %>%
        system
}


find_fastq_n_reads <- function(fastq){
    fastq %>% 
        str_c("wc -l ", .) %>% 
        system(intern = T) %>% 
        str_split(" ") %>% 
        .[[1]] %>% 
        .[[1]] %>% 
        as.integer %>% 
        divide_by(4) %>% 
        as.integer
}

create_fastq_sample_commands <- function(df, other_args = ""){
    pmap_chr(
        list(
            df$prefix,
            df$n_samples,
            df$p1_fastq_file,
            df$p2_fastq_file),
        create_fastq_sample_command,
        other_args)
}

create_fastq_sample_command <- function(
    prefix, n_samples, fastq_file1, fastq_file2 = "", other_args = ""){
    
    str_c(
        "fastq-sample",
        "-n", as.character(n_samples),
        "-o", prefix,
        "-r",
        other_args,
        fastq_file1,
        fastq_file2,
        sep = " ")
}

create_fastq_merge_command <- function(files, prefix){
    files %>% 
        str_c(collapse = " ") %>% 
        str_c("cat ", ., " > ", prefix, ".fastq")
}