### bioconductor install

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos='http://mirrors.ustc.edu.cn/CRAN/')
BiocManager::install()

### install packages

library(BiocManager)

install(c('ggtree','tidytree','treeio','Biostrings','phylotools','DECIPHER'))


