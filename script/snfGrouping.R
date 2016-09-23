###########################################
# Similarity Network Fusion analysis of HGSC datasets for subtyping
# 
# Chang, T
# ~~~~~~~~~~~~~
# Output subtype grouping for 2-8 subtypes
# Dataset input required

args <- commandArgs(trailingOnly = T)
# args <- c(2, 8, 20, 0.5, 20, "TCGA_HumanMethylation27", "TCGA_miRNA_HiSeq")

############################
# 0. Parse Arguments
############################
CLUSTERMIN <- as.numeric(args[1])
CLUSTERMAX <- as.numeric(args[2])
K_VALUE <- as.numeric(args[3])
SIGMA <- as.numeric(args[4])
T_VALUE <- as.numeric(args[5])
DATA1 <- args[6]
DATA2 <- args[7]

############################
# 1. Load Libraries
############################
library(SNFtool)

############################
# 2. Load Data
############################
Exp1 <- read.table(paste("data/raw/", DATA1, sep=""), header = T, stringsAsFactors = F, row.names = 1)
Exp2 <- read.table(paste("data/raw/", DATA2, sep=""), header = T, stringsAsFactors = F, row.names = 1)

############################
# 3. Sample Validation and Removal
############################
# Sample comparison between datasets
Exp1Col <- colnames(Exp1)
Exp2Col <- colnames(Exp2)
ExpColFilter <- Exp1Col[Exp1Col %in% Exp2Col]

# Filter samples and NA values in dataset 1
Exp1Filter <- Exp1[,ExpColFilter]
Exp1Filter.has.na <- apply(Exp1Filter, 1, function(x){any(is.na(x))})
Exp1Filter <- Exp1Filter[!Exp1Filter.has.na,]

# Filter samples and NA values in dataset 2
Exp2Filter <- Exp2[,ExpColFilter]
Exp2Filter.has.na <- apply(Exp2Filter, 1, function(x){any(is.na(x))})
Exp2Filter <- Exp2Filter[!Exp2Filter.has.na,]

############################
# 4. SNF Grouping
############################
SNFData <- list(Exp1Filter = Exp1Filter, Exp2Filter = Exp2Filter)
SNFGroup <- list()
WFused <- list()

for (i in CLUSTERMIN:CLUSTERMAX) {
  WData <- list()
  
  for (j in 1:length(SNFData)) {
    
    # Calculate pairwise distances
    distanceMatrix <- dist2(as.matrix(t(SNFData[[j]])), as.matrix(t(SNFData[[j]])))
    # Construct similarity graphs
    WData[[j]] <- affinityMatrix(distanceMatrix, K_VALUE, SIGMA)
    # Test1 <- affinityMatrix(distanceMatrix, K_VALUE, SIGMA)
  }
  
  # Fuse all graphs
  WFused[[i]] = SNF(WData, K=K_VALUE, t=T_VALUE)
  # Test2 = SNF(WData, K=K_VALUE, t=T_VALUE)
  
  # Perform spectral clustering
  SNFGroup[[i]] <- spectralClustering(WFused[[i]], i)
  # Test3 <- spectralClustering(Test2, 2)
}

############################
# 5. Data Output
############################
# Convert groupings into a labeled data frame
SNFGroupMatrix <- data.frame(matrix(unlist(SNFGroup), nrow=length(SNFGroup[[CLUSTERMIN]]), byrow=F))
rownames(SNFGroupMatrix) <- ExpColFilter
SNFColName <- c()
for (i in 1:length(SNFGroupMatrix)) {
  SNFColName[i] <- paste("K=", CLUSTERMIN + (i-1), sep = "")
}
colnames(SNFGroupMatrix) <- SNFColName

# Write grouping data to CSV
fName <- paste("SNF_",DATA1,"-",DATA2,".csv",sep="")
write.csv(SNFGroupMatrix, file = file.path("output","groups",fName), row.names = TRUE)

# Write filtered datasets to CSV
fName <- paste(DATA1, ".csv", sep = "")
write.csv(Exp1Filter, file = file.path("data","filtered",fName), row.names = TRUE)

fName <- paste(DATA2, ".csv", sep = "")
write.csv(Exp2Filter, file = file.path("data","filtered",fName), row.names = TRUE)