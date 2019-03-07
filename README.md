# DataCodes_iOmicsPASS

The data and codes are used in the following paper:

H.K.W. Koh, D. Fermin, C. Vogel, K.P. Choi, R. Ewing, and H. Choi, iOmicsPASS: network-based integration of multi-omics data for predictive subnetwork discovery. bioRxiv 374520; doi: https://doi.org/10.110/374520 (2018)

## Application: Breast Invasive Carcinoma (TCGA)
This folder contains 2 sub-folders, one for the analysis using 3 types of -omics data (CNV, mRNA and Protein) and the other for the analysis using 2 type of data (mRNA and Protein).

### Sub-folder: Analysis using CNV
This folder contains the input parameter file as well as the R codes ("Rcodes" folder) used to generate the heatmaps and plots in the manuscript in the analysis where all three types of -omics data (i.e. CNV, mRNA and Protein) were used.

### Sub-folder: Analysis without CNV
This folder contains the input parameter file as well as the R codes ("Rcodes" folder) used to generate the heatmaps and plots in the manuscript in the analysis where only two types of -omics data (i.e. mRNA and Protein) were used.

## Network Data
This folder contains two types of biological network files, (1) Transcription factor regulatory network and (2) Protein-protein interaction network, that were manually curated from various expermentally validated sources (see manuscript for more details). The two files were used as input in iOmicsPASS to link the genes with protein level information in the application of invasive breast cancer - TCGA. The folder also contain a (3) Pathway module information that is curated from consensusPath database and Gene Ontology (GO) and used as input in iOmicsPASS for subnetwork enrichment module.


## Datasets
Due to the large file size, we were unable to upload some of the files onto github. Hence, here we share the link to the dropbox folder containing all the necessary files to run the R-codes appended:
https://www.dropbox.com/sh/4yemvkk0lgz4y17/AACqnpiL0LOrtVNO73lCw-LUa?dl=0

The folder holds processed data files of each -omics (i.e. after removing low quality genes and data imputation), the group information file for the samples across all the datasets and as the output files of running iOmicsPASS separatetly for both analysis with or without CNV.


