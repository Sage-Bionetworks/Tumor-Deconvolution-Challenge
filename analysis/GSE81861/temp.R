library(tidyverse)
library(synapseClient)
library(data.table)
library(magrittr)
library(MCPcounter)
library(gplots)
library(heatmap.plus)
library(RColorBrewer)
library(ggfortify)

tmp_dir  <- "/home/aelamb/tmp/"

annotation_id <- "syn11898281"
count_id       <- "syn11867912"

setwd(tmp_dir)
synapseCacheDir(tmp_dir)
synapseLogin()

create_df_from_synapse_id <- function(syn_id, downloadLocation = NULL, unzip = F){
    path <- download_from_synapse(syn_id, downloadLocation)
    if(unzip) path <- str_c("zcat ", path)
    path %>% 
        fread %>% 
        as_data_frame 
}

download_from_synapse <- function(syn_id, downloadLocation = NULL){
    path = synGet(syn_id, downloadLocation = downloadLocation)@filePath
    return(path)
}

format_gene_strings <- function(strings){
    c(strings) %>%
        str_match("^([:print:]+):([:digit:]+)-([:digit:]+)_([:print:]+)_(ENSG[R]*[:digit:]+).[0-9]+$") %>% 
        .[,-1] %>% 
        split(1:nrow(.)) %>% 
        map_chr(str_c, collapse = ";;;")
}

format_column_names <- function(strings){
    strings %>% 
        str_split("__") %>% 
        map_chr(function(lst) lst[[1]])
}

transpose_df <- function(df, var1, var2){
    df %>% column_to_rownames(var1) %>% 
        as.matrix %>% 
        t %>% 
        data.frame %>% 
        rownames_to_column(var2)
}

get_count_summary_by_group <- function(samples, count_matrix, fn){
    m <- count_matrix[,samples]
    if(length(samples) == 1) return(m)
    apply(m, 1, fn)
}

annotation_df <- create_df_from_synapse_id(annotation_id) %>% 
    select(title, characteristics_ch1) %>% 
    separate(characteristics_ch1, c("string", "person"), ": ") %>% 
    mutate(person = str_replace_all(person, ";", "")) %>% 
    select(title, person) %>% 
    set_names(c("sample", "person"))


count_df <- count_id %>% 
    create_df_from_synapse_id("./", T) %>% 
    mutate(parsed_column = format_gene_strings(V1)) %>% 
    separate(parsed_column, into = c("chr", "start", "end", "Hugo", "Ensembl"), sep = ";;;", remove = T) %>% 
    mutate(gene_length = abs(as.numeric(start) - as.numeric(end)) -1) %>% 
    select(-c(V1, chr, start, end)) %>% 
    select(Hugo, Ensembl, gene_length, everything())

gene_metadata_df <- select(count_df, Hugo, Ensembl, gene_length) 


gene_metadata_df$Hugo[[which(gene_metadata_df$Ensembl == "ENSG00000261040")]] <- "WFDC21P" 


sample_metadata_df <- count_df %>% 
    names %>% 
    .[-c(1:3)] %>% 
    data_frame("sample_id_string" = .) %>% 
    separate(sample_id_string, into = c("sample", "cell_type"), sep = "__") %>% 
    left_join(annotation_df) %>% 
    mutate(combined_name = str_c(person, "_", cell_type))

color_df1 <- data_frame(
    "cell_type" = unique(sample_metadata_df$cell_type),
    "color_cell_type" = brewer.pal(nlevels(as.factor(sample_metadata_df$cell_type)), name = 'Dark2'))

color_df2 <- data_frame(
    "person" = unique(sample_metadata_df$person),
    "color_person" = brewer.pal(nlevels(as.factor(sample_metadata_df$person)), name = 'Set3'))

sample_metadata_df <- sample_metadata_df %>% 
    left_join(color_df1) %>% 
    left_join(color_df2)
    


count_df2 <- count_df
names(count_df)[-c(1:3)] <- sample_metadata_df$sample

genes = read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)$`HUGO symbols`


count_matrix <- count_df %>% 
    select(-c(gene_length, Hugo)) %>% 
    column_to_rownames("Ensembl") %>% 
    as.matrix

mean_count_m <- sample_metadata_df %>%
    split(.$combined_name) %>%
    map(use_series, sample) %>%
    map(get_count_summary_by_group, count_matrix, mean) %>%
    do.call("cbind", .) %>%
    set_rownames(gene_metadata_df$Hugo)

mean_rpm_m <- apply(mean_count_m, 2, function(x) (x/sum(x))*1000000) %>% 
    .[rownames(.) %in% genes,]

person_metadata_df <- sample_metadata_df %>% 
    select(-sample) %>% 
    distinct %>% 
    arrange(person, cell_type)

color_matrix <- person_metadata_df %>% 
    select(color_cell_type, color_person) %>% 
    set_colnames(c("Cell_Type", "Person")) %>% 
    as.matrix
    
    
heatmap.plus(log2(mean_rpm_m + 1), 
          ColSideColors = color_matrix)

heatmap.plus(mean_rpm_m, 
             ColSideColors = color_matrix)

x <- t(mean_rpm_m)
autoplot(prcomp(x), data = person_metadata_df, colour = 'person')
autoplot(prcomp(x), data = person_metadata_df, colour = 'cell_type')


# count_matrix2 <- count_df2 %>% 
#     filter(Hugo %in% genes) %>% 
#     select(-c(gene_length, Ensembl)) %>% 
#     column_to_rownames("Hugo") %>% 
#     as.matrix %>% 
#     add(1) %>% 
#     log2

# 
# 
# log_count_df <- count_matrix2 %>% 
#     t %>% 
#     data.frame 
#     
# 

# 
# write.table(count_matrix2, "./counts.tsv")



# x <- d3heatmap(count_matrix2)
# 
# svg("./heatmap.svg", width = 100, height = 100)
# heatmap.2(count_matrix2)
# dev.off()

# 

# 
# median_count_m <- annotation_df %>% 
#     filter(sample %in% sample_metadata_df$sample) %>% 
#     split(.$person) %>% 
#     map(use_series, sample) %>% 
#     map(get_count_summary_by_person, count_matrix, median) %>% 
#     do.call("cbind", .) %>% 
#     set_rownames(gene_metadata_df$Hugo) 
# 
# mean_count_m <- annotation_df %>% 
#     filter(sample %in% sample_metadata_df$sample) %>% 
#     split(.$person) %>% 
#     map(use_series, sample) %>% 
#     map(get_count_summary_by_person, count_matrix, mean) %>% 
#     do.call("cbind", .) %>% 
#     set_rownames(gene_metadata_df$Hugo) 
# 
# 
# 
# median_rpm_m <- apply(median_count_m, 2, function(x) (x/sum(x))*1000000)
# median_rpkm_m <- (median_rpm_m / gene_metadata_df$gene_length) * 1000
# median_rpk_m <- median_count_m / gene_metadata_df$gene_length 
# median_tpm_m <- apply(median_rpk_m, 2, function(x) (x/sum(x))*1000000)
# 
# 
# mean_rpm_m <- apply(mean_count_m, 2, function(x) (x/sum(x))*1000000)
# mean_rpkm_m <- (mean_rpm_m / gene_metadata_df$gene_length) * 1000

# mean_rpk_m <- mean_count_m / gene_metadata_df$gene_length
# mean_tpm_m <- apply(mean_rpk_m, 2, function(x) (x/sum(x))*1000000)
# 
# result <- MCPcounter.estimate(median_rpm_m, featuresType = "HUGO_symbols")
# result2 <- MCPcounter.estimate(median_rpkm_m, featuresType = "HUGO_symbols")
# result3 <- MCPcounter.estimate(median_rpm_m, featuresType = "HUGO_symbols")
# result4 <- MCPcounter.estimate(median_tpm_m, featuresType = "HUGO_symbols")
# 
# 
# result5 <- MCPcounter.estimate(mean_rpm_m, featuresType = "HUGO_symbols")
# result6 <- MCPcounter.estimate(mean_rpkm_m, featuresType = "HUGO_symbols")
# result7 <- MCPcounter.estimate(mean_rpm_m, featuresType = "HUGO_symbols")
# result8 <- MCPcounter.estimate(mean_tpm_m, featuresType = "HUGO_symbols")
# 
# res_list <- list(result, result2, result3, result4, result5, result6, result7, result8)
# res_names <- c("median_rpm", "median_rpkm", "median_rpk", "median_tpm", "mean_rpm", "mean_rpkm", "mean_rpk", "mean_tpm")
# 
# 
# format_result_matrix <- function(matrix){
#     matrix %>% 
#         round(2) %>% 
#         t %>% 
#         data.frame %>% 
#         rownames_to_column("sample") %>% 
#         gather(cell_type, p, T.cells:Fibroblasts)
# }
# 
# 
# result_df <- res_list %>% 
#     map(format_result_matrix) %>% 
#     map2(res_names, function(df, name) set_names(df, c("person", "cell_type", name))) %>% 
#     reduce(left_join) %>% 
#     arrange(person)
#     
# 
# cell_type_p_df <- left_join(sample_metadata_df, annotation_df) %>% 
#     group_by(person, cell_type) %>% 
#     summarise(count = n()) %>%
#     mutate(p = (count / sum(count) * 100)) %>% 
#     select(-count)
# 
# 
# 
# 
# median_count_m2 <-  median_count_m[rownames(median_count_m) %in% genes$`HUGO symbols`,]
# 
# 
# ExampleEstimates=MCPcounter.estimate(MCPcounterExampleData,featuresType="affy133P2_probesets")
# 
