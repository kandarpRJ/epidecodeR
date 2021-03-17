# epidecodeR

epidecodeR, an R package capable of integrating DNA/RNA epigenetic data generated from a host of epigenomic or epitranscriptomic techniques such as ChIP-seq, ATAC-seq, m6A-seq, etc. and dysregulated gene lists in the form of differential gene expression, ribosome occupancy or differential protein translation and identify impact of dysregulation of genes caused due to varying degrees of DNA/RNA chemical modifications associated with the genes. epidecodeR generates cumulative distribution function (CDF) plots showing shifts in trend of overall log2FC between genes divided into groups based on the degree of modification associated with the genes. The tool also tests for significance of difference in log2FC between groups of genes.

Installation procedure  #It is recommended to run following steps in Rstudio after installing R and Rstudio

1) Install R

Choose the CRAN mirror nearest to you (https://cran.r-project.org/mirrors.html)

2) Install Rstudio # Recommended but not necessary 

Install Rstudio for your OS (https://rstudio.com/products/rstudio/)

3) Install bioconductor and devtools packages required for installing epidecodeR and its dependencies # epidecodeR depends on following packages, which must be installed prior to installing epidecodeR. NOTE: except rtracklayer and GenomicRanges other packages will install automatically during installation of epidecodeR.

Dependencies:<br/>
  a. EnvStats<br/>
  b. ggplot2<br/>
  c. rtracklayer<br/>
  d. GenomicRanges<br/>
  e. rstatix<br/>
  f. ggpubr<br/>

if (!require("BiocManager")) {<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;install.packages("BiocManager")<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;library(BiocManager)<br/>
}
<br/>
if (!require("rtracklayer")){<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BiocManager::install("rtracklayer")<br/>
}
<br/>
if (!require("GenomicRanges")){<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BiocManager::install("GenomicRanges")<br/>
}
<br/>
if (!require("devtools")) {<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;install.packages("devtools")<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;library(devtools)<br/>
}
<br/>

4) Install epidecodeR from github (https://github.com/kandarpRJ/epidecodeR)

install_github("kandarpRJ/epidecodeR")<br/>
library (epidecodeR)
