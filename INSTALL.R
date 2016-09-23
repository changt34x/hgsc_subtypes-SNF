###########################################
# Similarity Network Fusion analysis of HGSC datasets for subtyping
# 
# Chang, T
# ~~~~~~~~~~~~~
# Install required packages

mirror <- "http://cran.us.r-project.org"

# Required CRAN packages
cran_packages <- C(
  'sfsmisc',
  'cluster',
  'ggplot2',
  'gplots',
  'RColorBrewer',
  'reshape'
)

install.packages(cran_packages, repos = mirror)

# Required Bioconductor packages
source("https://bioconductor.org/biocLite.R")
bioc_packages <- c(
  'Biobase',
  'SNFTool'
)

biocLite(bioc_packages, suppressUpdates = TRUE)