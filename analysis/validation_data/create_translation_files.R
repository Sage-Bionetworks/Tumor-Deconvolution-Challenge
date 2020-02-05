library(biomaRt)
synapser::synLogin()
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
bm <- getBM(
    attributes = c(
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "hgnc_symbol"
    ), 
    mart = ensembl
)
readr::write_tsv(bm, "translation.tsv")  
synapser::File("translation.tsv", parent = "syn21571479") %>% 
    synapser::synStore()
file.remove("translation.tsv")
