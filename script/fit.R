###########################################
# Similarity Network Fusion analysis of HGSC datasets for subtyping
# 
# Chang, T
# ~~~~~~~~~~~~~
# Compute subtype grouping statistics
# 
# References:
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B, Konecny, G., Goode, E., C.S., Doherty, J.A. (2016).
# Unpublished: Cross-population analysis of high-grade serious ovarian cancer does not support four subtypes.

args <- commandArgs(trailingOnly = T)
# args <- c(2, 8, 20, 20, "TCGA_HumanMethylation27", "TCGA_miRNA_HiSeq")

############################
# 0. Parse Arguments
############################
CLUSTERMIN <- as.numeric(args[1])
CLUSTERMAX <- as.numeric(args[2])
ITERMAX <- as.numeric(args[3])
STARTS <- as.numeric(args[4])
DATA1 <- args[5]
DATA2 <- args[6]

############################
# 1. Load Libraries
############################
library(sfsmisc)
library(cluster)
library(ggplot2)

# Load functions
source("script/functions/functions.R")

############################
# 2. Load Data
############################
Exp1 <- read.csv(paste("data/filtered/", DATA1, ".csv", sep=""), header = T, stringsAsFactors = F, row.names = 1)
Exp2 <- read.csv(paste("data/filtered/", DATA2, ".csv", sep=""), header = T, stringsAsFactors = F, row.names = 1)

ExpData <- list(Exp1, Exp2)

############################
# 3. Silhouette Plots
############################
Name <- paste("SNF_",DATA1,"-",DATA2,".csv",sep="")

RunSilhouette(ExpData[[2]], Name, 2, 4)
