# LUAD

This repository contains matadata and codes necessary for analysis of the scRNA data of human lung adenocarcinoma (LUAD) presented in Han and Sinjab etal. Nature, in revision (2022). 

Cells with low complexity libraries or likely cell debris have been removed.
Batch effect was corrected by Harmony when analyzing the no-malignant cells. 

Please contact hkadara@mdanderson.org and LWang22@mdanderson.org with any questions.

## Downloading the data

Cell level metadata as well as the further level quantification data is available in the /Input_data folder, with detailed annotation in the corresponding codes parts

Sequencing data for P1-P5 were previously generated and deposited in the European Genome–phenome Archive (EGA) under the accession number EGAS00001005021. Sequencing data generated in this study including human scRNA-seq data (P6 – P16) and human spatial transcriptomics data will be deposited in EGA under the same accession number (EGAS00001005021) or be deposited in NCBI GEO (mouse scRNA-seq data and mouse spatial transcriptomics data) and made publicly available with the final version of this article.

### Requirements

##Tested on macOS Monterey (12.6.1)

##R packages:
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

-InferCNV(v1.3.2)
-CytoTRACE(v0.3.3)
