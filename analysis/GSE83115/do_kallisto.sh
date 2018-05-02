git clone https://github.com/Sage-Bionetworks/kallisto_cwl
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
cwltool kallisto_cwl/index.cwl --index_file_string GRCH38.idx --fasta_file Homo_sapiens.GRCh38.cdna.all.fa.gz &> index_log.txt &
rm Homo_sapiens.GRCh38.cdna.all.fa.gz
Rscript create_kallisto_quants.R &> kallisto_log.txt &