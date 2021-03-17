# epidecodeR

epidecodeR, an R package capable of integrating DNA/RNA epigenetic data generated from a host of epigenomic or epitranscriptomic techniques such as ChIP-seq, ATAC-seq, m6A-seq, etc. and dysregulated gene lists in the form of differential gene expression, ribosome occupancy or differential protein translation and identify impact of dysregulation of genes caused due to varying degrees of DNA/RNA chemical modifications associated with the genes. epidecodeR generates cumulative distribution function (CDF) plots showing shifts in trend of overall log2FC between genes divided into groups based on the degree of modification associated with the genes. The tool also tests for significance of difference in log2FC between groups of genes.

Installation procedure  #It is recommended to run following steps in Rstudio after installing R and Rstudio

1) Install R

Choose the CRAN mirror nearest to you (https://cran.r-project.org/mirrors.html)

2) Install Rstudio # Recommended but not necessary 

Install Rstudio for your OS (https://rstudio.com/products/rstudio/)

3) Install bioconductor and devtools packages required for installing epidecodeR and its dependencies # epidecodeR depends on following packages, which must be installed prior to installing epidecodeR. NOTE: except rtracklayer and GenomicRanges other packages will install automatically during installation of epidecodeR.

Dependencies:
a. EnvStats
b. ggplot2
c. rtracklayer
d. GenomicRanges
e. rstatix
f. ggpubr

if(!require("BiocManager")) {<br/>
&nbsp&nbspinstall.packages(c("BiocManager", "devtools"))<br/>
&nbsp&nbsplibrary(BiocManager)<br/>
}
<br/>
if(!require("rtracklayer")){<br/>
&nbsp&nbspBiocManager::install(c("rtracklayer", "GenomicRanges"))<br/>
}
<br/>
4) Install epidecodeR from github (https://github.com/kandarpRJ/epidecodeR)

library(devtools)
install_github("kandarpRJ/epidecodeR")
library (epidecodeR)
