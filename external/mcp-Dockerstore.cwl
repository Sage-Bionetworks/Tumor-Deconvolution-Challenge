#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [Rscript, /usr/local/bin/mcp_counter.R]

doc: "run MCPcounter"

requirements:
- class: InlineJavascriptRequirement

hints:
  DockerRequirement:
    dockerPull: quay.io/cri-iatlas/iatlas-tool-mcpcounter 

inputs:

  input_expression_file:
    type: File
    inputBinding:
      position: 1
      prefix: --input_expression_file
    doc: Path to input matrix of microarray expression data. Tab separated file with features in rows and samples in columns
    
  output_file_string:
    type: string
    inputBinding:
      prefix: --output_file
    default: "./output_file.tsv"
    doc: path to write output file

  features_type:
    type: string
    inputBinding:
      prefix: --features_type
    default: "affy133P2_probesets"
    doc: Type of identifiers for expression features. Defaults to 'affy133P2_probesets' for Affymetrix Human Genome 133 Plus 2.0 probesets. Other options are 'HUGO_symbols' (Official gene symbols) or 'ENTREZ_ID' (Entrez Gene ID)

  input_probeset_file:
    type: ["null", File]
    inputBinding:
      prefix: --input_probeset_file
    doc: Path to input table of gene data. Tab separated file of probesets transcriptomic markers and corresponding cell populations. Fetched from github by a call to read.table by default, but can also be a data.frame

  input_gene_file:
    type: ["null", File]
    inputBinding:
      prefix: --input_gene_file
    doc: Path to input table of gene data. Tab separated file of genes transcriptomic markers (HUGO symbols or ENTREZ_ID) and corresponding cell populations. Fetched from github by a call to read.table by default, but can also be a data.frame

outputs:

  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_file_string)
    doc: see output_string



