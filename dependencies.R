# first, install Bioconductor for Cells application
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
# next, install main dependencies
deps <- c("foreach", "doParallel", "dbscan", "ggplot2", "mclust","xtable", "plyr",
          "dplyr", "reshape2", "grid", "latex2exp", "gridExtra", "plyr", "dplyr",
          #"scran", "scuttle", "M3Drop",
          "latex2exp", "cowplot", "uwot",
          "RColorBrewer", "RSSL", "mvtnorm",
          "ggforce", "KODAMA", "devtools")
install.packages(deps)
rm(deps)
# install Bioconductor packages
BiocManager::install("scran")
BiocManager::install("scuttle")
BiocManager::install("M3Drop")
# install mcclust.ext from Github
devtools::install_github("sarawade/mcclust.ext")
# install fold from Github
devtools::install_github("adombowsky/FOLD/fold")
