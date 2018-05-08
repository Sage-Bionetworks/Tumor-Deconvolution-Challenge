git clone https://github.com/Sage-Bionetworks/kallisto_cwl
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
cwltool kallisto_cwl/index.cwl --index_file_string GRCH38.idx --fasta_file Homo_sapiens.GRCh38.cdna.all.fa.gz &> index_log.txt &
rm Homo_sapiens.GRCh38.cdna.all.fa.gz
mkdir tmp
Rscript Tumor-Deconvolution-Challenge/analysis/GSE83115/dl_fastqs_for_kallisto.R.R &> dl_log.txt &
Rscript Tumor-Deconvolution-Challenge/analysis/GSE83115/create_kallisto_quants.R &> kallisto_log.txt &
Rscript Tumor-Deconvolution-Challenge/analysis/GSE83115/ul_kallisto_quants.R &> ul_log.txt &