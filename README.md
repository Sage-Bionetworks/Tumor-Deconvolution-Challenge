# Tumor-Deconvolution-Challenge

See analysis/validation-analysis/run-deconvolution-method-on-challenge-data.R, which describes how to run a deconvolution method against the Challenge data, using xCell as an example, and to compare it to the Challenge results.

The following scripts in the github repository 
https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge
were used to generate the indicated tables and figures:

| Table number | Table name | Script that produced the table |
| --- | --- | --- |
| S1 | analysis/validation-analysis/figs/rerun-validation-mixture-and-distribution-effects.tsv | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| S2 | analysis/sample-level-analysis/figs/figs/sample-level-comparison.tsv | analysis/sample-level-analysis/plot-sample-level-summaries.R |
| S10 | training-sample-metadata/geo-expression-array-immune-cells.csv | training-sample-metadata/download-and-format-training-sample-metadata.R | 
| S11 | training-sample-metadata/geo-rnaseq-immune-cells.csv | training-sample-metadata/download-and-format-training-sample-metadata.R |
| 1 | method-annotation-table/deconv-method-description-table-1.pdf | method-annotation-table/deconv-method-description-table-1.tex |
| 2 | method-annotation-table/deconv-method-description-table-2.pdf | method-annotation-table/deconv-method-description-table-2.tex |
| S12 | analysis/timing/method-run-times.tsv | analysis/timing/format-timing-results.R |
| S13 | training-sample-metadata/SampleDetails.kallisto.csv.gz | training-sample-metadata/download-and-format-training-sample-metadata.R |
| S14 | training-sample-metadata/DA505_training_metadata.csv | https://figshare.com/s/944b97b01c91763e8dcd |
| S15 | training-sample-metadata/detailedMetadataOfSamples.xlsx | training-sample-metadata/download-and-format-training-sample-metadata.R |
| S16 | training-sample-metadata/Metadata.csv | training-sample-metadata/download-and-format-training-sample-metadata.R |
| S19 | analysis/cancer-validation/figs/cancer-validation-correlations.tsv | analysis/cancer-validation/score-cancer-datasets.R |
| S20 | analysis/cancer-validation/figs/cancer-validation-dataset-comparison-pvals.tsv | analysis/cancer-validation/score-cancer-datasets.R |

| Figure number | File name | Script that produced the figure |
| --- | --- | --- |
| Fig S1 | analysis/summary/purified-samples-marker-heatmap-protein-coding-genes-no-title.png | analysis/summary/plot-marker-heatmap.R |
| Fig 2 |  analysis/validation-analysis/figs/fig-validation-round-1-performance.png | analysis/validation-analysis/perform-bootstrap-rerun-validation.R, plot-bootstrap-rerun-validation.R |
| Fig S2 | external-figs/da_505-supp-fig1.png | DA505 writeup: https://www.synapse.org/#!Synapse:syn20674744/wiki/603943 |
| Fig S3 | external-figs/Biogem_ComparisonDataset.png | Biogem writeup: https://www.synapse.org/#!Synapse:syn20551146/wiki/594139 |
| Fig S4 | analysis/validation-analysis/figs/fig-validation-performance-across-rounds.png | make-validation-performance-figs.R |
| Fig S5 | analysis/validation-analysis/figs/fig-validation-all-performance.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig 3 | analysis/validation-analysis/figs/fig-validation-round-1-strip-and-heatmap-merged-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S6 | analysis/validation-analysis/figs/fig-validation-round-1-merged-strip-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S7  | analysis/validation-analysis/figs/fig-validation-heatmap-round-1-coarse-and-fine-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S8 | analysis/validation-analysis/figs/fig-validation-round-1-coarse-and-fine-strip-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S9 | analysis/validation-analysis/figs/fig-validation-heatmap-rounds-2-and-3-merged-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S10 | analysis/validation-analysis/figs/fig-validation-round-2-merged-strip-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S11 | analysis/validation-analysis/figs/fig-validation-heatmap-round-2-coarse-and-fine-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S12 | analysis/validation-analysis/figs/fig-validation-round-2-coarse-and-fine-strip-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S13 | analysis/validation-analysis/figs/fig-validation-round-3-merged-strip-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S14 | analysis/validation-analysis/figs/fig-validation-heatmap-round-3-coarse-and-fine-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig S15 | analysis/validation-analysis/figs/fig-validation-round-3-coarse-and-fine-strip-cell-type.png | analysis/validation-analysis/plot-bootstrap-rerun-validation.R |
| Fig 4 | analysis/sample-level-analysis/figs/sample-level-metric-swarm-round-1.png | analysis/sample-level-analysis/plot-sample-level-summaries.R |
| Fig 5 | analysis/specificity-analysis/figs/spillover-summary.png | analysis/specificity-analysis/plot-validation-spillover-results.R |
| Fig S16 | analysis/specificity-analysis/figs/spillover-all-scores-coarse-grained-round-1.png | analysis/specificity-analysis/plot-validation-spillover-results.R |
| Fig S17 | analysis/specificity-analysis/figs/spillover-all-scores-fine-grained-round-1.png | analysis/specificity-analysis/plot-validation-spillover-results.R |
| Fig 6 | analysis/in-silico-admixtures/figs/sensitivity-spikein-and-summary.png | analysis/in-silico-admixtures/plot-sensitivity-results.R |
| Fig S18 | analysis/cancer-validation/figs/fig-cancer-validation-heatmap-wu-and-pelka.png | analysis/cancer-validation/score-cancer-datasets.R |
| Fig S19 | analysis/cancer-validation/figs/fig-cancer-validation-heatmap-all.png | analysis/cancer-validation/score-cancer-datasets.R |
| Fig 7 | analysis/cancer-validation/figs/fig-cancer-validation-per-cell-type.png | analysis/cancer-validation/score-cancer-datasets.R |
| Fig S20 | analysis/cancer-validation/figs/fig-cancer-validation-legend-edited.png | analysis/cancer-validation/score-cancer-datasets.R |
| Fig S21 | analysis/summary-table/figs/binned-score-heatmap-score-sorted-ignore-nas-comparator.png | analysis/summary-table/make-summary-table.R |



