# epigenome methylation profiling by array Illumina Infinium 
# Yao, S., Zhao, L., Chen, S., Wang, H., Gao, Y., Shao, N. Y., Dai, M., & Cai, H. (2023). Cervical cancer immune infiltration microenvironment identification, construction of immune scores, assisting patient prognosis and immunotherapy. 
# Frontiers in immunology, 14, 1135657. https://doi.org/10.3389/fimmu.2023.1135657

GSE30759_matrix <- getGEO("GSE30759", GSEMatrix = TRUE)
GSE30759_phenoData <- pData(GSE30759_matrix[[1]])
