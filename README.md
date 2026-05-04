# Blood-based proteomic signature of brain amyloidosis

# Table of contents
* [Introduction](#introduction)
* [Content](#content)
* [Data](#data)
* [Requirements](#requirements)
* [License](#license)
* [Instructions](#instructions)

# Introduction
This repository contains the code for bioinformatics analyses described in the article "Blood-based proteomic signature of amyloidosis: Identification of novel regulators of amyloid load".

This project investigated plasma proteomics data from the SomaScan v4.1 assay to identify proteins associated with cerebral amyloidosis and link them to clinical variability within Aβ+ individuals. Idnetified proteins were leveraged to create a weigted proteomic score that correlates with amyloid burden, classical AD biomarkers, and clinical disease severity. Proteomic clustering and biological pathway analyses were performed to understand underlying biology.

# Content
The code covers the following main analysis steps:

1. Data pre-processing
2. Differential abundance analysis (DAA)
3. Proteomic score calculation, trait association, and clustering
4. Pathway enrichment analyses
   
# Data
Proteomics data analysed in this study is available at:
- Knight-ADRC: https://live-knightadrc-washu.pantheonsite.io/professionals-clinicians/request-center-resources/
- Bio-Hermes: Members of the global research community can place a data use request via the AD Discovery Portal (https://discover.alzheimersdata.org/). Access is contingent upon adherence to the Bio-Hermes Data Use Agreement.

# Requirements
The code was written in R (version 4.3.0) and relies on multiple R and Bioconductor packages, including:
- dplyr (1.1.4)
- tidyr (1.3.1)
- data.table (1.16.4)
- pheatmap (1.0.12)
- clusterProfiler (4.12.6)
- ReactomePA (1.48.0)
- DOSE (3.30.5)
- enrichplot (1.24.4)
- ggplot2 (3.5.2)
- EnhancedVolcano (1.22.0)

- Additional packages listed at the beginning of each R script

# License
The code is available under the MIT License.

# Instructions
The code was tested on R 4.3.0 on Linux operating systems, but should be compatible with later versions of R installed on current Linux, Mac, or Windows systems.

To run the code, the correct working directory containing the input data must be specified at the beginning of the R-scripts, otherwise the scripts can be run as-is.

The scripts should be run in the following order:

    data_preparation.R

    differential_expression_analysis.R

    trait_association_clustering.R

    proteomic_score_and_clinical_association.R

    pathway_enrichment_analysis.R

