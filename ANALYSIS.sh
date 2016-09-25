#!/bin/bash

exec &>hgsc_analysis.out 

############################################
# Similarity Network Fusion analysis of HGSC datasets for subtyping
# 
# Chang, T
# ~~~~~~~~~~~~~
# This script stores instructions to perform subtyping of HGSC samples. 
# All scripts and relevant files are included in the repository and the 
# workflow depends on the running sequential scripts within  the larger
# folder structure. See the README for general information and INSTALL.R
# for package dependencies.
# ~~~~~~~~~~~~~~~~~~~~~
############################################

################
# Part 0:
# Constants
################
CLUSTERMIN=2
CLUSTERMAX=4
DATA1="TCGA_HumanMethylation27
DATA2="TCGA_miRNA_HiSeq
K=20
SIGMA=0.5
T=20
ITERMAX=20
STARTS=20


################
# Part 1:
# SNF Grouping
################
R --no-save --args $CLUSTERMIN $CLUSTERMAX $K $SIGMA $T $DATA1 $DATA2 \
< script/snfGrouping.R

################
# Part 2:
# Fit Testing
################
R --no-save --args $CLUSTERMIN $CLUSTERMAX $ITERMAX $STARTS $DATA1 $DATA2 \
< script/fit.R

################
# Part 3:
# Correlation Matrix
################
R --no-save --args $CLUSTERMIN $CLUSTERMAX $DATA1 $DATA2 \
< script/correlationMatrix.R