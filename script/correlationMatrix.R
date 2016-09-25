###########################################
# Similarity Network Fusion analysis of HGSC datasets for subtyping
# 
# Chang, T
# ~~~~~~~~~~~~~
# Create correlation matrix heatmaps
# 
# References:
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B, Konecny, G., Goode, E., C.S., Doherty, J.A. (2016).
# Unpublished: Cross-population analysis of high-grade serious ovarian cancer does not support four subtypes.

args <- commandArgs(trailingOnly = T)
# args <- c(2, 8, "TCGA_HumanMethylation27", "TCGA_miRNA_HiSeq")

############################
# 0. Parse Arguments
############################
CLUSTERMIN <- as.numeric(args[1])
CLUSTERMAX <- as.numeric(args[2])
DATA1 <- args[3]
DATA2 <- args[4]

krange <- CLUSTERMIN:CLUSTERMAX

############################
# 1. Load Libraries
############################
library(gplots)
library(RColorBrewer)
library(reshape)

# Load heatmap.3 function
source("script/functions/heatmap3.R")

# Load functions
source("script/functions/functions.R")

############################
# 2. Load Data
############################
Exp1 <- read.csv(paste("data/filtered/", DATA1, ".csv", sep=""), header = T, stringsAsFactors = F, row.names = 1)
Exp2 <- read.csv(paste("data/filtered/", DATA2, ".csv", sep=""), header = T, stringsAsFactors = F, row.names = 1)

ExpData <- list(Exp1, Exp2)

############################
# 3. Build Correlation Matrices
############################
name <- paste("SNF_",DATA1,"-",DATA2,sep="")
membfile <- list.files(path = "output/groups", pattern = paste(Name, ".csv", sep=""))
ClusterAssign <- read.csv(file = file.path("output", "groups", membfile), row.names = 1)

for (centroid in 1:length(krange)) {
  dataclustertmp <- CorMatrixOrder(ExpData[[1]], ClusterAssign, as.numeric(paste(centroid)))
  
  colorPool <- c('skyblue1', 'tomato')
  if (centroid == 2) {
    colorPool <- c(colorPool, "springgreen")
  } else if (centroid == 3) {
    colorPool <- c(colorPool, "springgreen", "violet")
  }
  
  if (name == "GSE26712_eset" & centroid == 1) {
    colorPool <- c("gray75", "gray60")
  } else if (name == "GSE26712_eset" & centroid == 2) {
    colorPool <- c("gray75", "gray60", "gray45")
  } else if (name == "GSE26712_eset" & centroid == 3) {
    colorPool <- c("gray75", "gray60", "gray45", "gray30")
  }
  
  sideColors <- as.matrix(rep(colorPool, as.vector(table(ClusterAssign[ ,centroid]))))
  
  heatFile <- paste(DATA1, krange[centroid], sep = "_")
  heatTitle <- c("Correlation Matrix\n", paste(unlist(strsplit(name, "_"))[1], "k =", 
                                               krange[centroid], sep = " "))
  
  color_range <- colorRampPalette(c("blue", "White", "red"))(n = 2999)
  
  mid_cut = mean(dataclustertmp) - (3 * sd(dataclustertmp))
  mid_cut_high = quantile(dataclustertmp, probs = .65)
  col_breaks = c(seq(min(dataclustertmp), mid_cut, length = 100), seq(mid_cut,mid_cut_high, length = 1600), 
                 seq(mid_cut_high, 1, length = 1300))
  
  png(file.path("output","figures",paste(heatFile, ".png", sep="")), width = 1000, height = 900)
  
  heatmap.3(dataclustertmp, symm = T, trace = 'none', Rowv = NULL, Colv = "Rowv", dendrogram = 'none', 
            key = F, labRow = F, labCol = F, col = color_range, breaks = col_breaks, main = heatTitle, 
            ColSideColors = sideColors, ColSideColorsSize = 2, RowSideColors = t(sideColors), 
            RowSideColorsSize = 2)
  
  dev.off()
}

for (centroid in 1:length(krange)) {
  dataclustertmp <- CorMatrixOrder(ExpData[[2]], ClusterAssign, as.numeric(paste(centroid)))
  
  colorPool <- c('skyblue1', 'tomato')
  if (centroid == 2) {
    colorPool <- c(colorPool, "springgreen")
  } else if (centroid == 3) {
    colorPool <- c(colorPool, "springgreen", "violet")
  }
  
  if (name == "GSE26712_eset" & centroid == 1) {
    colorPool <- c("gray75", "gray60")
  } else if (name == "GSE26712_eset" & centroid == 2) {
    colorPool <- c("gray75", "gray60", "gray45")
  } else if (name == "GSE26712_eset" & centroid == 3) {
    colorPool <- c("gray75", "gray60", "gray45", "gray30")
  }
  
  sideColors <- as.matrix(rep(colorPool, as.vector(table(ClusterAssign[ ,centroid]))))
  
  heatFile <- paste(DATA2, krange[centroid], sep = "_")
  heatTitle <- c("Correlation Matrix\n", paste(unlist(strsplit(name, "_"))[1], "k =", 
                                               krange[centroid], sep = " "))
  
  color_range <- colorRampPalette(c("blue", "White", "red"))(n = 2999)
  
  mid_cut = mean(dataclustertmp) - (3 * sd(dataclustertmp))
  mid_cut_high = quantile(dataclustertmp, probs = .65)
  col_breaks = c(seq(min(dataclustertmp), mid_cut, length = 100), seq(mid_cut,mid_cut_high, length = 1600), 
                 seq(mid_cut_high, 1, length = 1300))
  
  png(file.path("output","figures",paste(heatFile, ".png", sep="")), width = 1000, height = 900)
  
  heatmap.3(dataclustertmp, symm = T, trace = 'none', Rowv = NULL, Colv = "Rowv", dendrogram = 'none', 
            key = F, labRow = F, labCol = F, col = color_range, breaks = col_breaks, main = heatTitle, 
            ColSideColors = sideColors, ColSideColorsSize = 2, RowSideColors = t(sideColors), 
            RowSideColorsSize = 2)
  
  dev.off()
}