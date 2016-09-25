#####################################
# FUNCTIONS
#####################################

silhouette_SimilarityMatrix<-function(group, similarity_matrix)
{
  similarity_matrix=as.matrix(similarity_matrix)
  similarity_matrix<-(similarity_matrix+t(similarity_matrix))/2
  diag(similarity_matrix)=0
  normalize <- function(X) X / rowSums(X)
  similarity_matrix<-normalize(similarity_matrix)
  
  n <- length(group)
  if(!all(group == round(group))) stop("'integer")
  cluster_id <- sort(unique(group <- as.integer(group)))
  k <- length(cluster_id)
  if(k <= 1 || k >= n)
    return(NA)
  doRecode <- (any(cluster_id < 1) || any(cluster_id > k))
  if(doRecode)
    group <- as.integer(fgroup <- factor(group))
  cluster_id <- sort(unique(group))
  
  wds <- matrix(NA, n,3, dimnames =list(names(group), c("cluster","neighbor","sil_width")))  
  for(j in 1:k)
  { 
    index <- (group == cluster_id[j])
    Nj <- sum(index)
    wds[index, "cluster"] <- cluster_id[j]
    dindex <- rbind(apply(similarity_matrix[!index, index, drop = FALSE], 2,
                          function(r) tapply(r, group[!index], mean)))
    maxC <- apply(dindex, 2, which.max)
    wds[index,"neighbor"] <- cluster_id[-j][maxC]
    s.i <- if(Nj > 1) {
      a.i <- colSums(similarity_matrix[index, index])/(Nj - 1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      ifelse(a.i != b.i, (a.i - b.i) / pmax(b.i, a.i), 0)
    } else 0
    wds[index,"sil_width"] <- s.i
  }
  attr(wds, "Ordered") <- FALSE
  class(wds) <- "silhouette"
  wds
}

# Determine specific silhouette widths for each dataset for each k=2 to k=4
RunSilhouette <- function(dataset, Name, incr) {
  # ~~~~~~~~~~~~~~
  # Output silhouette width plots
  #
  # Args: 
  # dataset: A gene expression matrix with genes as rows and samples as columns
  # Name: The name of the eset 
  # 
  # Returns:
  # Silhouette width plots
  # 
  # References:
  # Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B, Konecny, G., Goode, E., C.S., Doherty, J.A. (2016).
  # Unpublished: Cross-population analysis of high-grade serious ovarian cancer does not support four subtypes.
  # ~~~~~~~~~~~~~~
  
  dataExprs <- data.frame(t(dataset))
  
  dataDist <- dist(dataExprs)
  
  # Read in subtype groupings
  membfile <- list.files(path = "output/groups", pattern = paste(Name, ".csv", sep=""))
  ClusterAssign <- read.csv(file = file.path("output", "groups", membfile), row.names = 1)
  
  # Perform silhouette width analyses
  K2_sil <- silhouette(as.numeric(paste(ClusterAssign$K.2)), dataDist)
  K3_sil <- silhouette(as.numeric(paste(ClusterAssign$K.3)), dataDist)
  K4_sil <- silhouette(as.numeric(paste(ClusterAssign$K.4)), dataDist)
  
  # Output images
  png(file.path("output","figures",paste(Name, incr, "_Silhouette.png", sep="")), width = 1100, height = 800)
  
  # Get appropriate margins
  par(mfrow = c(1, 3))
  par(mar = c(2, 1.5, 1, 1.5))
  par(cex.axis = 0.7, cex = 2, cex.main = 2)
  
  plot(K2_sil, main = "", xlab = "", sub = "", col = c('skyblue1', 'tomato'), do.col.sort = T, do.n.k = F)
  plot(K3_sil, main = "", xlab = "", sub = "", col = c('skyblue1', 'tomato', 'springgreen'), do.col.sort = T, do.n.k = F)
  plot(K4_sil, main = "",  xlab = "", sub = "", col = c('skyblue1', 'tomato', 'springgreen', 'violet'), do.col.sort = T, do.n.k = F)
  dev.off()
}

CorMatrixOrder <- function (DataSet, ClusterMemb, ClusterColumn) {  
  # ~~~~~~~~~~~~~~
  # Outputs a dataframe of cluster membership for k = 3 and k = 4 using Global MAD genes
  #
  # Args: 
  # DataSet: A gene expression matrix with genes as rows and samples as columns
  # ClusterMemb: a dataframe with cluster assignments
  # ClusterColumn: which column the cluster memberships are stored
  # 
  # Returns:
  # The order by which the samples in each correlation matrix heatmap are plotted
  # ~~~~~~~~~~~~~~
  
  # Initialize an order for the correlation matrix
  DataSet_order <- c()
  
  # Determine the number of unique clusters in the given cluster column
  uniqueClusters <- length(unique(ClusterMemb[ ,ClusterColumn]))
  
  # For each cluster, perform hierarchical clustering to determine the order, 
  # within each cluster assignment, for presentation in the correlation matrix
  for (clus in 1:uniqueClusters) {
    # Select the gene expression of the samples with the cluster assignment
    cluster.subset <- DataSet[ ,ClusterMemb[ ,ClusterColumn] == clus]
    
    # Observe the distance matrix for these samples
    cluster.dist <- dist(t(cluster.subset))
    
    # Perform heirarchical clustering on this distance matrix
    cluster.hclust <- hclust(cluster.dist)
    
    # Order the subset of samples based on the order resulting from the heirarchical clustering
    cluster.order <- cluster.subset[ ,cluster.hclust$order]
    cluster.order <- as.matrix(cluster.order)
    
    # Reorder the gene expression matrices based on first the k means clusters and next based on 
    # the heirarchical clusters
    DataSet_order <- cbind(DataSet_order, cluster.order)
  }
  
  # Obtain the correlation matrix and return
  corReady <- cor(DataSet_order)
  return(corReady)
}