# High-Grade Serous Ovarian Cancer Subtypes - Similarity Network Fusion Analysis

## Summary
This repository contains the process used to perform Similarity Network Fusion analysis
for subtyping on high-grade serous ovarian cancer samples. The sample data used for the
procedure  was The Cancer Genome Atlas' Ovarian Cancer data set, obtainable from Xena,
but the procedure is designed to be adaptable for data from a variety of studies. 

## Analysis
For best results, the analysis should be performed on  Ubuntu 16.04 with R 3.3.1. 

Prior to analysis, run INSTALL.R to obtain the newest packages required.

Place datasets (e.g. on miRNA Expression and Methylation) in /data/raw.

Update the file names in ANALYSIS.sh accordingly.

## Acknowledgements
This work was highly influenced by "High-Grade Serious Ovarian Cancer Subtypes - Why has
the field settled on four?" available [here](https://github.com/greenelab/hgsc_subtypes).