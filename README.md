# PING
 Probabilistic inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data

This is R implementation of statitical method proposed in my paper

Zhang X, Robertson G, Woo S, Hoffman B, Gottardo R (2012). “Probabilistic Inference for Nucleosome Positioning with MNase-Based or Sonicated Short-Read Data.” PLoS ONE, 7. DOI:10.1371/journal.pone.0032095.

Please cite this paper, if you used this R package in your research. Thanks!

To install this R package, please use the following R code:

    library(devtools)
    install_github("ubcxzhang/PING")


This is a clone of my Bioconductor package PING version 2.28.0 https://www.bioconductor.org/packages/release/bioc/html/PING.html

So, this package can be installed from Bioconductor using  the following R code:

    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("PING")
