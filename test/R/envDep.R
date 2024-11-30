### Bioconductor download
install.packages("BiocManager")

### Install required packages
BiocManager::install(c("GSEABase","GSVA","msigdbr","clusterProfiler"), ask = F, update = F)
BiocManager::install(c("GEOquery","limma","impute"), ask = F, update = F)
BiocManager::install(c("org.Hs.eg.db","org.Mm.eg.db"), ask = F, update = F)
BiocManager::install(c("DESeq2","edgeR"), ask = F, update = F)
BiocManager::install("enrichplot", ask = F, update = F)
BiocManager::install("devtools", ask = F, update = F)
BiocManager::install("WGCNA", ask = F, update = F)
BiocManager::install("data.table", ask = F, update = F)
BiocManager::install("tximport", ask = F, update = F)
BiocManager::install("tidyverse", ask = F, update = F)
BiocManager::install("DOSE", ask = F, update = F)
BiocManager::install("patchwork", ask = F, update = F)
BiocManager::install("RBGL", ask = F, update = F) # Dependency for Vennerable
BiocManager::install("pathview", ask = F, update = F)
BiocManager::install(c("STRINGdb","ggraph","igraph"), ask = F, update = F)
install.packages("Vennerable", repos="http://R-Forge.R-project.org") # Install Vennerable package
install.packages("statmod") # Install some other basic packages
options()$repos
install.packages(c("FactoMineR", "factoextra"))
install.packages(c("ggplot2", "pheatmap","ggpubr","ggthemes",
                   "ggstatsplot","ggsci","ggsignif"))
install.packages("rvcheck")
(.packages())  # View currently loaded packages

# Update all packages
rvcheck::update_all(check_R = F, which = c("CRAN", "BioC", "github"))
