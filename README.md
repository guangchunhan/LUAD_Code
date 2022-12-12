# LUAD

This repository contains matadata and codes necessary for analysis of the scRNA data of human lung adenocarcinoma (LUAD) presented in Han and Sinjab etal. Nature, in revision (2022). 

Cells with low complexity libraries or likely cell debris have been removed.
Batch effect was corrected by Harmony when analyzing the no-malignant cells. 

Please contact hkadara@mdanderson.org and LWang22@mdanderson.org with any questions.

## Downloading the data

Cell level metadata is available in the provided /Input_data/Cell metaData.rda, which contains QC statistics of cells, clinical information of associated samples and cell types. 

All the raw data are uploading to EGA at the moment and the accession number will be released very soon (xxx-xx-xx).

## Data visualization

### Requirements

Tested on macOS Big Sur

1. R version: 4.1.2
2. R packages
   - ggplot2
   - data.table
   - Seurat
   - dplyr
   - tidyr
   - ggpubr
   - RColorBrewer
   - monocle
   - smoother
   - pheatmap
   - Hmisc
   - monocle3
   - ggrepel

3. inferCNV_version(xxx)
4. CytoTRACE_version(xxx)
