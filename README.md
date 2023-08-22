# Two Sample MR workflow for Terra.bio

Import the workflow to your Terra workspace using the link below.

- [Dockstore](https://dockstore.org/workflows/github.com/anand-imcm/terra-TwoSampleMR-wf1:main?tab=info)

Locate the 'Launch with' widget at the top right of the Dockstore workflow page, and select the 'Terra' platform option. 

## About

A WDL-based workflow that utilizes the R package "TwoSampleMR" to generate single SNP and multi SNP forest plots. This workflow is designed to facilitate Mendelian randomization analysis using [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/) package and visualize the results in the form of forest plots.

## WDL tasks

- `get_all_exposures` [__Task 1__] The initial task in the WDL workflow extracts a list of exposures from the input exposure summary statistics file. The list is divided into user-defined chunks for further processing.
- `twosamplemr`  [__Task 2__] The second task in the WDL script utilizes the `scatter` function to execute an R script in parallel for each exposure set obtained from Task 1.
  - `src/two_sample_mr__script_1.R` [Rscript] A script which an exposure summary statistics file, an outcome summary statistics file, and a list of exposures to test.
- `combine_objects` [__Task 3__] The third task in the workflow collects the output files generated in Task 2 and runs an R script as mention below.
  - `src/gather_mr_outputs__script_2.R` [Rscript] The script aggregates all the objects generated in the previous step into several `.rds` files.
- `generate_plot` [__Task 4__]  The final step that generates the forest plots.
  - `src/nv_wf1_multiANDsingle_v1.r` [Rscript] The custom Rscript is designed to generate two different plots, one for multi-SNPs and another for single-SNPs , and save them in the PDF format.

## Data preparation

In this step, the necessary input data for the analysis is prepared. This may include genetic variant data, exposure and outcome data, and any other required information. The input data should be formatted appropriately for compatibility with the main Rscript.

- Exposure summary statistics [File] A tab-delimited file with eight required columns: SNP, beta, se, effect_allele, other_allele, eaf, pval, Marker. __Marker__ column is required. Example: `data/proteomics_summary_data.FINAL.ALL.txt`
- Outcome summary statistics [File] A tab-delimited file with seven required columns: SNP, beta, se, effect_allele, other_allele, eaf, pval. Example: `data/Kunkle_etal_2019_IGAP_Summary_statistics.with_allelefreqs.FINAL.ALL.txt`

## Inputs

- Exposure summary statistics [tsv]
- Outcome summary statistics [tsv]
- clumping [Float] nClumping filter. Default: `0.01`
- npval [Float] nPval exposure filter. Default: `'5e-08'`
- chunk_size [Int] The parallel processing of markers will occur in subsets of user defined chunk sizes. Default: `50`
- fdr [Float] Threshold for the fdr adjusted p values used for the analyses. Default: `0.05`
- nsplit [Int] The maximum number of exposures shown in each page of the pdf output file. Default: `10`

## Output

- single SNP joined plot [pdf]
- multi SNP joined plot [pdf]
- rds objects for downstream analysis

## Components

- Docker images
  - debian:stable-20230502-slim
  - mrcieu/twosamplemr:0.5.6
  - ghcr.io/anand-imcm/vis:v1.0
- R packages
  - littler
  - optparse
  - qpdf
  - berryFunctions
  - ggforestplot
