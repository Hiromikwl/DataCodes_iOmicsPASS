# DataCodes_iOmicsPASS

The data and codes are used in the following paper:

H.K.W. Koh, D. Fermin, C. Vogel, K.P. Choi, R. Ewing, and H. Choi, iOmicsPASS: network-based integration of multi-omics data for predictive subnetwork discovery. bioRxiv 374520; doi: https://doi.org/10.110/374520 (2018)

## Application: Breast Invasive Carcinoma (TCGA)
This folder contains 3 sub-folders, one for the analysis using 3 types of -omics data (CNV, mRNA and Protein), one for the analysis using 2 type of data (mRNA and Protein) and the third folder containing the processed input files used for running iOmicsPASS.

### Sub-folder: Analysis using CNV
This folder contains the input parameter file used in iOmicsPASS using all three types of -omics data (i.e. CNV, mRNA and Protein) in the analysis. It also contains the direct output files ("Output" folder), as well as the R codes ("Rcodes" folder) used to generate the heatmaps and plots in the manuscript.

### Sub-folder: Analysis without CNV
This folder contains the input parameter file used in iOmicsPASS using only two types of -omics data (i.e. mRNA and Protein) in the analysis. It also contains the direct output files ("Output" folder), as well as the R codes ("Rcodes" folder) used to generate the corresponding heatmaps and plots in the manuscript.

### Sub-folder: InputFiles_BRCA
This folder contains the processed data files of each -omics, after removing low quality genes and data imputation, as well as the group information for the samples across all the datasets.


## Network Data
This folder contains two types of biological network files, (1) Transcription factor regulatory network and (2) Protein-protein interaction network, that were manually curated from various expermentally validated sources (see manuscript for more details). The two files were used as input in iOmicsPASS to link the genes with protein level information in the application of invasive breast cancer - TCGA. The folder also contain a (3) Pathway module information that is curated from consensusPath database and Gene Ontology (GO) and used as input in iOmicsPASS for subnetwork enrichment module.
