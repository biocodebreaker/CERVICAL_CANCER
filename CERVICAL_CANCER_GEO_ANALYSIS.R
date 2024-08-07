BiocManager::install("GEOquery")

library(GEOquery)

library("dplyr")

library("tidyverse")


# Nie, H., Bu, F., Xu, J., Li, T., & Huang, J. (2020). 29 immune-related genes pairs signature predict the prognosis of cervical cancer patients. 
# Scientific reports, 10(1), 14152. https://doi.org/10.1038/s41598-020-70500-5

# gse44001 sample size = 300

gse44001_matrix <- getGEO("gse44001", GSEMatrix = TRUE)
gse44001_phenoData <- pData(gse44001_matrix[[1]])

gse44001_phenoData_data_processing <- gse44001_phenoData[[22]]
# The software package to extract raw data was Genomestudio. 
# The intensity of the probes was transformed by binary logarithm 
# and then was normalized using the quantile normalization method

gse44001_phenoData_platform_id <- gse44001_phenoData[[23]]

# GPL14951 Illumina HumanHT-12 WG-DASL V4.0 R2 expression beadchip


#Get Annotation/ Feature data
gse44001_matrix_featureData <- fData(gse44001_matrix[[1]])



# GSE39001

GSE39001_matrix <- getGEO("GSE39001", GSEMatrix = TRUE)

GSE39001_matrix_phenoData <- pData(GSE39001_matrix[[1]])


GSE39001_matrix_phenoData_subset <- pData(GSE39001_matrix[[1]]) %>%
  as.data.frame() %>%
  select(sample_GSM = geo_accession,
         sampleType = `clinic:ch1`,
         sampleID = `sample id:ch1`,
         tissueType = `tissue:ch1`)

GSE63514_matrix_phenoData_subset <- pData(GSE63514_matrix[[1]]) %>%
  as.data.frame() %>%
  select(sample_GSM = geo_accession,
         sampleType = `title`,
         tissueType = `tissue type:ch1`)




# Microarrays were normalized using the RMA algorithm (Robust Multichip Average)
# in the Affymetrix Expression Console.

GSE39001_matrix_phenoData_data_processing <- GSE39001_matrix_phenoData[[20]]

GSE39001_matrix_phenoData_platform_id <- GSE39001_matrix_phenoData[[21]]

# GPL201 [HG-Focus] Affymetrix Human HG-Focus Target Array



# GSE63514

GSE63514_matrix <- getGEO("GSE63514", GSEMatrix = TRUE)

GSE63514_phenoData <- pData(GSE63514_matrix[[1]])

# GeneChip RMA (GC-RMA) is an improved form of RMA that is able to use 
# the sequence-specific probe affinities of the GeneChip probes to attain more accurate gene expression values.
GSE63514_phenoData_data_processing <- GSE63514_phenoData[[23]]

GSE63514_phenoData_platform_id <- GSE63514_phenoData[[24]]
#  GPL570 Affymetrix Human Genome U133 Plus 2.0 Array



#Get ExpressionSet data from Series Matrix file (GSE)

GSE39001_expression_data <- exprs(GSE39001_matrix[[1]])

GSE63514_expression_data <- exprs(GSE63514_matrix[[1]])

#First transform the expression matrix to a data.frame 
# library("tidyverse")

exprs(GSE39001_matrix[[1]]) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") -> GSE39001_expression_data

dim(GSE39001_expression_data)
# [1] 8793   56

exprs(GSE63514_matrix[[1]]) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") -> GSE63514_expression_data

dim(GSE63514_expression_data)
# [1] 54675   129

# Get Annotation/ Feature data
GSE39001_featureData <- fData(GSE39001_matrix[[1]])

GSE63514_featureData <- fData(GSE63514_matrix[[1]])


# Join the Feature Data with the Expression data

# subset the GSE39001_featureData data frame to include only the columns "ID", "ENTREZ_GENE_ID" and "Gene Symbol" 
GSE39001_featureData_subset <- GSE39001_featureData %>%
  select("ID", "Gene Symbol", "ENTREZ_GENE_ID")


GSE39001_expression_data %>%
  inner_join(GSE39001_featureData_subset) -> GSE39001_expression_data_GeneSymbol_original

# subset the GSE63514_featureData  data frame to include only the columns "ID" and "Gene Symbol" 
GSE63514_featureData_subset <- GSE63514_featureData %>%
  select("ID", "Gene Symbol")

# Please note that inner_join `by = join_by(ID)` so Gene Symbol column appears last

GSE63514_expression_data %>%
  inner_join(GSE63514_featureData_subset) -> GSE63514_expression_data_GeneSymbol_original


# Rearrange the columns so that Gene Symbol column comes first
GSE39001_expression_data_GeneSymbol_original_arranged <- GSE39001_expression_data_GeneSymbol_original %>%
  select("Gene Symbol", ENTREZ_GENE_ID, everything())

GSE63514_expression_data_GeneSymbol_original_arranged <- GSE63514_expression_data_GeneSymbol_original %>%
  select("Gene Symbol", everything())


# First remove the ENTREZ_GENE_ID column
GSE39001_expression_data_GeneSymbol_original_filtered <- GSE39001_expression_data_GeneSymbol_original_arranged %>%
  select(-ENTREZ_GENE_ID)

# Rename the "Gene Symbol" column to "GeneSymbol"

GSE39001_expression_data_GeneSymbol_original_sorted <- GSE39001_expression_data_GeneSymbol_original_sorted %>%
  rename(GeneSymbol = "Gene Symbol")
dim(GSE39001_expression_data_GeneSymbol_original_sorted)
# [1] 8793   57

GSE63514_expression_data_GeneSymbol_original_arranged <- GSE63514_expression_data_GeneSymbol_original_arranged %>%
  rename(GeneSymbol = "Gene Symbol")

# First remove the ID column
GSE63514_expression_data_GeneSymbol_original_arranged <- GSE63514_expression_data_GeneSymbol_original_arranged %>%
  select(-ID)

# Sort the data frame by the "Gene Symbol" column

GSE39001_expression_data_GeneSymbol_original_sorted <- GSE39001_expression_data_GeneSymbol_original_filtered %>%
  arrange(`Gene Symbol`)

GSE63514_expression_data_GeneSymbol_original_arranged_sorted <- GSE63514_expression_data_GeneSymbol_original_arranged  %>%
  arrange(`GeneSymbol`)

# Remove rows with blank or NA values in GeneSymbol column
GSE39001_expression_data_GeneSymbol_original_sorted_filtered <- GSE39001_expression_data_GeneSymbol_original_sorted %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "")

dim(GSE39001_expression_data_GeneSymbol_original_sorted_filtered)
# [1] 8688   57

# Remove rows with blank or NA values in GeneSymbol column
GSE63514_expression_data_GeneSymbol_original_arranged_sorted_filtered <- GSE63514_expression_data_GeneSymbol_original_arranged_sorted  %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "")

dim(GSE63514_expression_data_GeneSymbol_original_arranged_sorted_filtered)
# [1] 45782   129


# First remove the ID column
GSE39001_expression_data_GeneSymbol_sorted_filtered <- GSE39001_expression_data_GeneSymbol_original_sorted_filtered %>%
  select(-ID)

dim(GSE39001_expression_data_GeneSymbol_sorted_filtered)
# [1] 8688   56

# Calculate Median Expression per Gene Symbol for GSE39001

# Remove the first column (GeneSymbol) for aggregation
df_for_aggregation <- GSE39001_expression_data_GeneSymbol_sorted_filtered[, -1]
df_for_aggregation <- GSE63514_expression_data_GeneSymbol_original_arranged_sorted_filtered[, -1]


# Calculate the median for each gene symbol

gene_symbol_medians <- aggregate(. ~ GSE39001_expression_data_GeneSymbol_sorted_filtered$GeneSymbol, data = df_for_aggregation, median)
gene_symbol_medians <- aggregate(. ~ GSE63514_expression_data_GeneSymbol_original_arranged_sorted_filtered$GeneSymbol, data = df_for_aggregation, median)


# Rename the columns for clarity
colnames(gene_symbol_medians) <- c("GeneSymbol", colnames(df_for_aggregation))

GSE63514_Median_Expression_23520_129 <- gene_symbol_medians

dim(GSE39001_Median_Expression_limma_input)
#[1] 8524   56

GeneCards_MSigDB_Oxidative_Stress_genes_9570 <- readRDS("GeneCards_MSigDB_Oxidative_Stress_genes_9570.rds")

selected_genesymbols <- GeneCards_MSigDB_Oxidative_Stress_genes_9570$GeneSymbol

length(selected_genesymbols)
# [1] 9569

common_genes <- intersect(GSE39001_OS_genes_Median_Expression_limma_input$GeneSymbol, GeneCards_MSigDB_Oxidative_Stress_genes_9570$GeneSymbol)
length(common_genes)
# [1] 4439

# Filter rows in GSE63514_Median_Expression_23520_129 where "GeneSymbol" is in selected_genesymbols

GSE63514_OS_genes_Median_Expression_limma_input_7006 <- GSE63514_Median_Expression_23520_129[GSE63514_Median_Expression_23520_129$GeneSymbol %in% selected_genesymbols, ]
dim(GSE39001_OS_genes_Median_Expression_limma_input_4439)
# [1] 4439   56


# Convert data.table to data.frame
# TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal <- as.data.frame(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal)

class(GSE39001_OS_genes_Median_Expression_limma_input_4439)
# [1] "data.frame"

# set rownames to GeneSymbol
rownames(GSE39001_OS_genes_Median_Expression_limma_input_4439) <- GSE39001_OS_genes_Median_Expression_limma_input_4439$GeneSymbol
GSE39001_OS_genes_Median_Expression_limma_input_4439$GeneSymbol <- NULL

GSE39001_OS_genes_Median_Expression_limma_input_4439_12_Normal_43_Tumor <- GSE39001_OS_genes_Median_Expression_limma_input_4439


rownames(GSE63514_OS_genes_Median_Expression_limma_input_7006) <- GSE63514_OS_genes_Median_Expression_limma_input_7006$GeneSymbol
GSE63514_OS_genes_Median_Expression_limma_input_7006$GeneSymbol <- NULL

GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor <- GSE63514_OS_genes_Median_Expression_limma_input_7006[, c(1:24, 101:128)]
 
dim(GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor)
# [1] 7006   52

library(limma)

# Create the design matrix for GSE39001_OS_genes_Median_Expression_limma_input_4439_12_Normal_43_Tumor

sample_groups <- factor(rep(c("Normal", "Tumor"), c(12, 43)))

# Create the design matrix for GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor

sample_groups <- factor(rep(c("Normal", "Tumor"), c(24, 28)))


# Convert the sample_groups to a factor
sample_groups <- factor(sample_groups, levels = c("Normal", "Tumor"))

# Create a data frame for the design matrix
design_matrix <- model.matrix(~0 + sample_groups)

# Set appropriate column names in the design matrix
colnames(design_matrix) <- c("Normal", "Tumor")

# Normalize using voom 

# v <- voom(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910, design_matrix)

# Fit the linear model: The use of eBayes or treat with trend=TRUE is known as the limma-trend method (Law et al, 2014; Phipson et al, 2016).
# recommended for microarray data (Gordon Smyth)

fit <- lmFit(GSE39001_OS_genes_Median_Expression_limma_input_4439_12_Normal_43_Tumor, design_matrix)

fit <- lmFit(GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor, design_matrix)


# Define the contrast matrix
contrast_matrix <- makeContrasts(Normal_vs_Tumor = Normal - Tumor, levels = design_matrix)

# Apply the contrast to the fit
fit_contrast <- contrasts.fit(fit, contrast_matrix)

# Perform empirical Bayes moderation of the standard errors
fit_eBayes <- eBayes(fit_contrast, trend = TRUE)


# Get the results of differential expression analysis
GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor_DEGs <- topTable(fit_eBayes, coef = "Normal_vs_Tumor", number = Inf, adjust.method = "fdr", sort.by = "P")

# TCGA_GTEx_v8_tpm_CESC_OS_DEGs_8253_fit_treat_robust <- topTreat(fit_treat_robust, coef = "Tumor_vs_Normal", number = Inf, adjust.method = "fdr", sort.by = "P")


write.table(GSE39001_OS_genes_4439_12_Normal_43_Tumor_DEGs, file = "GSE39001_OS_genes_4439_12_Normal_43_Tumor_DEGs.txt", sep = "\t", row.names = TRUE)

write.table(GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor_DEGs, file = "GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor_DEGs.txt", sep = "\t", row.names = TRUE)


GSE39001_OS_genes_4439_12_Normal_43_Tumor_DEGs <- readRDS("GSE39001_OS_genes_4439_12_Normal_43_Tumor_DEGs.rds")

TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910_DEGs <- readRDS("TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910_DEGs.rds")

GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor_DEGs <- readRDS("GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor_DEGs.rds")

library(limma)

# Get Significant DEGs for TCGA_GTEx adj.P.Val < 0.05
TCGA_GTEx_CESC_OS_306_Tumor_22_Normal_Significant_DEGs_5710 <- TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910_DEGs[TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910_DEGs$adj.P.Val < 0.05, ]

# Get Significant DEGs for GSE63514 adj.P.Val < 0.05

GSE63514_CESC_OS_24_Normal_28_Tumor_Significant_DEGs_2174 <- GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor_DEGs[GSE63514_OS_genes_limma_input_7006_24_Normal_28_Tumor_DEGs$adj.P.Val < 0.05, ]

# Get Significant DEGs for GSE39001 adj.P.Val < 0.05

GSE39001_CESC_OS_12_Normal_43_Tumor_Significant_DEGs_1643 <- GSE39001_OS_genes_4439_12_Normal_43_Tumor_DEGs[GSE39001_OS_genes_4439_12_Normal_43_Tumor_DEGs$adj.P.Val < 0.05, ]

# Find unique and intersection genes for TCGA-GTEx, GSE63514 and GSE39001

genes_df1 <- rownames(TCGA_GTEx_CESC_OS_306_Tumor_22_Normal_Significant_DEGs_5710)

genes_df2 <- rownames(GSE63514_CESC_OS_24_Normal_28_Tumor_Significant_DEGs_2174)

genes_df3 <- rownames(GSE39001_CESC_OS_12_Normal_43_Tumor_Significant_DEGs_1643)


genes_common_df1_df2_not_df3 <- genes_df1[genes_df1 %in% genes_df2 & !(genes_df1 %in% genes_df3)]
length(genes_common_df1_df2_not_df3)
 #[1] 1107

genes_common_df1_df3_not_df2 <- genes_df1[genes_df1 %in% genes_df3 & !(genes_df1 %in% genes_df2)]
length(genes_common_df1_df3_not_df2)
 #[1] 851

genes_common_df2_df3_not_df1 <- genes_df2[genes_df2 %in% genes_df3 & !(genes_df2 %in% genes_df1)]
length(genes_common_df2_df3_not_df1)
 # [1] 104
 
  # Genes common to all three data.frames genes_df1, genes_df2 and genes_df3
  common_genes <- genes_df1[genes_df1 %in% genes_df2 & genes_df1 %in% genes_df3]
  
  length(common_genes)
  # [1] 484
  
  # genes_df1 ONLY (TCGA_GTEx_CESC_OS_306_Tumor_22_Normal_Significant_DEGs_5710) 
  5710 - (1107+484+851)
  [1] 3268
 

  # genes_df2 ONLY (GSE63514_CESC_OS_24_Normal_28_Tumor_Significant_DEGs_2174)
  
  2174 - (104+484+1107)
  [1] 479
  
  # genes_df3 ONLY (GSE39001_CESC_OS_12_Normal_43_Tumor_Significant_DEGs_1643)
  
  1643 - (104+484+851)
  
  #[1] 204
  
  # Check logFC of the 484 common genes in TCGA_GTEx, GSE63514 and GSE39001 
  
  print(TCGA_GTEx_GSE63514_GSE39001_common_genes_484_logFC_sorted["ADAM9", ])
  #       logFC  AveExpr        t     P.Value   adj.P.Val         B
  # ADAM9 0.6055853 7.384453 3.123489 0.001945524 0.003219477 -3.165859
  
  print(GSE63514_CESC_OS_24_Normal_28_Tumor_Significant_DEGs_2174["ADAM9", ])
  
  print(GSE39001_CESC_OS_12_Normal_43_Tumor_Significant_DEGs_1643["ADAM9", ])
  
  # TCGA_GTEx_GSE63514_GSE39001_common_genes_484_logFC <- TCGA_GTEx_CESC_OS_306_Tumor_22_Normal_Significant_DEGs_5710[TCGA_GTEx_GSE63514_GSE39001_common_genes_484, ]
  
  # sort by rownames i.e. Gene Symbols 
  
  # TCGA_GTEx_GSE63514_GSE39001_common_genes_484_logFC_sorted <- TCGA_GTEx_GSE63514_GSE39001_common_genes_484_logFC[order(rownames(TCGA_GTEx_GSE63514_GSE39001_common_genes_484_logFC)), ]
  
  
  # Venn diagram showing for TCGA-GTEx, GSE63514 and GSE39001 
  
  
  # Install and load the VennDiagram package (if not already done)
  #install.packages("VennDiagram")
  library(VennDiagram)
  
  
  x <- list("Immune DEGs" = 1:565, "Stromal DEGs" = 426:1481)
  ggVennDiagram(x) + 
    scale_fill_gradient(low = "#0000FF", high = "#FFFF00")
  
  
  
  
  
  # Install and load the VennDiagram package if not already installed
  # install.packages("VennDiagram")
  library(VennDiagram)
  
  # Define the gene sets
  immune_genes <- c("CCR4", "FMO2", "SELP", "CYP1B1")
  stromal_genes <- c("CCR4", "FMO2", "SELP", "CYP1B1", "MAGEA11", "OMD", "CNR1", "CNKSR2")
  
  
  x <- list("Immune DEGs" = c("CCR4", "FMO2", "SELP", "CYP1B1") , "Stromal DEGs" = ("CCR4", "FMO2", "SELP", "CYP1B1", "MAGEA11", "OMD", "CNR1", "CNKSR2"))
  ggVennDiagram(x) + 
    scale_fill_gradient(low = "#0000FF", high = "#FFFF00")
  
  
  x <- list("Immune DEGs" = c("CCR4", "FMO2", "SELP", "CYP1B1"), 
            "Stromal DEGs" = c("CCR4", "FMO2", "SELP", "CYP1B1", "MAGEA11", "OMD", "CNR1", "CNKSR2"))
  
  # Plot the Venn diagram
  ggVennDiagram(x) + 
    scale_fill_gradient(low = "#0000FF", high = "#FFFF00")
  
  
  # Create a list with the gene sets
  gene_sets <- list("Immune DEGs" = immune_genes, "Stromal DEGs" = stromal_genes)
  
  # Plot the Venn diagram
  venn.plot <- venn.plot(gene_sets, category.names = c("Immune DEGs", "Stromal DEGs"))
  grid.draw(venn.plot)
  
  setwd("/cloud/home/r1816512/CERVICAL_CANCER")
  

  # use the t() function to transpose the data.frame after 
  # set the GeneSymbol column as row names
  
  rownames(TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors) <- TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors$GeneSymbol
  
  # subset the GeneSymbols such that the 484 common_genes remain 
  
  TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484 <- TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors[TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors$GeneSymbol %in% common_genes, ]
  
  # Remove the GeneSymbol column
  TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484$GeneSymbol <- NULL
  
  dim(TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484)
  # [1] 484 306
  
  # Transpose the data.frame
  TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484_transposed <- t(TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484)
  
  TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484_transposed <- as.data.frame(TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484_transposed)
  
  
  TCGA_CESC_survival_time_censored <- readRDS("TCGA_CESC_survival_time_censored.rds")
  
  TCGA_CESC_survival_time_censored <- as.data.frame(TCGA_CESC_survival_time_censored)
  
  rownames(TCGA_CESC_survival_time_censored) <- TCGA_CESC_survival_time_censored$submitter_id
  
  TCGA_CESC_survival_time_censored$submitter_id <- NULL
  
  
  TCGA_CESC_survival_time_censored$months <- TCGA_CESC_survival_time_censored$time / 30.44
  
  TCGA_CESC_survival_time_event_months <- TCGA_CESC_survival_time_censored[, c(4,5,6)]
  
  
  # keep only those rows where "time" > 0
  TCGA_CESC_survival_time_event_months_filtered_294 <- TCGA_CESC_survival_time_event_months[TCGA_CESC_survival_time_event_months$time > 0, ]
  
  TCGA_CESC_survival_time_event_months_filtered_294_sorted <- TCGA_CESC_survival_time_event_months_filtered_294[order(rownames(TCGA_CESC_survival_time_event_months_filtered_294)), ]
  
  
  # Specify the row names you want to print
  row_names_to_print <- c("TCGA-HM-A6W2-01A", "TCGA-HM-A6W2-06A", "TCGA-UC-A7PG-01A", "TCGA-UC-A7PG-06A")
  
  # Subset the data.frame based on the specified row names
  duplicate_rows <- TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484_transposed[row_names_to_print, ]
  
   duplicate_rows[,1:5]
   # ABCF1 ACAA1 ACAT1 ACOX2  ACP1
   # TCGA-HM-A6W2-01A  4294  1902   565    53  7660
   # TCGA-HM-A6W2-06A  6662  1772   846   256 16882
   # TCGA-UC-A7PG-01A  4706  2018  1127    53  1935
   # TCGA-UC-A7PG-06A  6979  1964  1174    92  2945
  
  # remove rows belonging to metastatic tumors 
  rows_to_remove <- c("TCGA-HM-A6W2-06A", "TCGA-UC-A7PG-06A")
  
  # remove rows belonging to Primary tumors 
  
  rows_to_remove <- c("TCGA-HM-A6W2-01A", "TCGA-UC-A7PG-01A")
  
  
  # Remove specified rows
  TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors <- TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484_transposed[setdiff(rownames(TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors_484_transposed), rows_to_remove), ]
  
  dim(TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors)
  # [1] 304 484
  
  # Get the row names
  row_names <- rownames(TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors)
  
  # Remove the part after the last hyphen
  edited_row_names <- sub("-[^-]+$", "", row_names)
  
  # Assign the edited row names back to the data frame
  rownames(TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors) <- edited_row_names
  

  # Add submitter_id as a column
  TCGA_CESC_survival_time_event_months_filtered_294_sorted$submitter_id <- rownames(TCGA_CESC_survival_time_event_months_filtered_294_sorted)
  TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors$submitter_id <- rownames(TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors)
  
  # rearrange such that submitter_id column comes first
  
  TCGA_CESC_survival_time_event_months_filtered_294_arranged <- TCGA_CESC_survival_time_event_months_filtered_294_sorted %>%
    select(submitter_id, everything())
  
  TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors_arranged <- TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors %>%
    select(submitter_id, everything())
  
  # Merge based on submitter_id
  merged_data <- TCGA_CESC_survival_time_event_months_filtered_294_arranged %>%
    inner_join(TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors_arranged, by = "submitter_id")
  
  

  
  
# Has Age, Gender, Tumor Stage I-IV
# 20 tumors, 8 Normal
GSE6791_matrix <- getGEO("GSE6791", GSEMatrix = TRUE)
GSE6791_phenoData <- pData(GSE6791_matrix[[1]])


GSE67522_matrix <- getGEO("GSE67522", GSEMatrix = TRUE)
GSE67522_matrix_phenoData <- pData(GSE67522_matrix[[1]])


# Zhao, K., Ma, Z., & Zhang, W. (2022). Comprehensive Analysis to Identify SPP1 as a Prognostic Biomarker in Cervical Cancer. 
# Frontiers in genetics, 12, 732822. https://doi.org/10.3389/fgene.2021.732822


GSE63514_matrix <- getGEO("GSE63514", GSEMatrix = TRUE)

GSE7803_matrix <- getGEO("GSE7803", GSEMatrix = TRUE)

GSE9750_matrix <- getGEO("GSE9750", GSEMatrix = TRUE)


# Get Clinical/phenotypic data

GSE63514_phenoData <- pData(GSE63514_matrix[[1]])

GSE7803_phenoData <- pData(GSE7803_matrix[[1]])

GSE9750_matrix_phenoData <- pData(GSE9750_matrix[[1]])

# Wu, B., & Xi, S. (2021). Bioinformatics analysis of differentially expressed genes and pathways in the development of cervical cancer. 
# BMC cancer, 21(1), 733. https://doi.org/10.1186/s12885-021-08412-4

GSE64217 <- getGEO("GSE64217")
GSE64217_phenoData <- pData(GSE64217[[1]])

# GSM1566487 <- getGEO("GSM1566487")


GSE138080 <- getGEO("GSE138080")

GSE138080_phenoData <- pData(GSE138080[[1]])

# GASTRIC CANCER ANALYSIS 

GSE84426_matrix <- getGEO("GSE84426", GSEMatrix = TRUE)
GSE84426_phenoData <- pData(GSE84426_matrix[[1]])




GSE84433_phenoData_subset <- pData(GSE84433_matrix[[1]]) %>%
  as.data.frame() %>%
  select(sample_GSM = geo_accession,
         age = `age:ch1`,
         event = `death:ch1`,
         time = `duration overall survival:ch1`,
         Nstage = `pnstage:ch1`,
         Tstage = `ptstage:ch1`,
         Gender = `Sex:ch1`)

GSE84426_phenoData_subset <- pData(GSE84426_matrix[[1]]) %>%
  as.data.frame() %>%
  select(sample_GSM = geo_accession,
         age = `age:ch1`,
         event = `death:ch1`,
         time = `duration overall survival:ch1`,
         Nstage = `pnstage:ch1`,
         Tstage = `ptstage:ch1`,
         Gender = `Sex:ch1`)

# dim(GSE84426_phenoData_subset)
# [1] 76  7


#Get ExpressionSet data from Series Matrix file (GSE)

#GSE84433_expression_data <- exprs(GSE84433[[1]])


#Join the Feature Data with the Expression data
#First transform the expression matrix to a data.frame 
#library("tidyverse")

exprs(GSE84433_matrix[[1]]) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") -> GSE84433_expression_data

exprs(GSE84426_matrix[[1]]) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") -> GSE84426_expression_data

# dim(GSE84426_expression_data)
# [1] 48708    77


#Get Annotation/ Feature data
GSE84433_featureData <- fData(GSE84433_matrix[[1]])

GSE84426_featureData <- fData(GSE84426_matrix[[1]])


#subset the GSE84433_featureData data frame to include only the columns "ID", "Entrez_Gene_ID" and "Symbol" 
GSE84433_featureData_subset <- GSE84433_featureData %>%
  select(ID, Entrez_Gene_ID, Symbol)


GSE84433_expression_data %>%
  inner_join(GSE84433_featureData_subset) -> GSE84433_expression_data_GeneSymbol_original

GSE84426_featureData_subset <- GSE84426_featureData %>%
  select(ID, Entrez_Gene_ID, Symbol)


GSE84426_expression_data %>%
  inner_join(GSE84426_featureData_subset) -> GSE84426_expression_data_GeneSymbol_original


#Check if genes to be validated are present 

genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")

genes_to_check %in% GSE84433_expression_data_GeneSymbol$Symbol
#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

genes_to_check %in% GSE84426_expression_data_GeneSymbol_original$Symbol

# Rearrange the columns
GSE84433_expression_data_GeneSymbol_arranged <- GSE84433_expression_data_GeneSymbol_original %>%
  select(Symbol, Entrez_Gene_ID, everything())

GSE84426_expression_data_GeneSymbol_arranged <- GSE84426_expression_data_GeneSymbol_original %>%
  select(Symbol, Entrez_Gene_ID, everything())

# Sort the data frame by the "Symbol" column
GSE84433_expression_data_GeneSymbol_sorted <- GSE84433_expression_data_GeneSymbol_arranged %>%
  arrange(Symbol)

GSE84426_expression_data_GeneSymbol_sorted <- GSE84426_expression_data_GeneSymbol_arranged %>%
  arrange(Symbol)

# Remove rows with blank or NA values in Symbol column
GSE84426_expression_data_GeneSymbol_filtered <- GSE84426_expression_data_GeneSymbol_sorted %>%
  filter(!is.na(Symbol) & Symbol != "")


# First remove the Entrez_Gene_ID column
GSE84433_expression_data_GeneSymbol_filtered <- GSE84433_expression_data_GeneSymbol_filtered %>%
  select(-Entrez_Gene_ID)

GSE84433_expression_data_GeneSymbol_filtered <- GSE84433_expression_data_GeneSymbol_filtered %>%
  select(-ID)

# Remove Entrez_Gene_ID and ILMN Probe ID
GSE84426_expression_data_GeneSymbol_filtered <- GSE84426_expression_data_GeneSymbol_filtered %>%
  select(-Entrez_Gene_ID, -ID)

#library(dplyr)

#library("tidyverse")

# Calculate Average Expression per Gene Symbol using the aggregate function

# Remove the first column (GeneSymbol) for aggregation
df_for_aggregation <- GSE84433_expression_data_GeneSymbol_filtered[, -1]

# Calculate the mean for each gene symbol
gene_symbol_means <- aggregate(. ~ GSE84433_expression_data_GeneSymbol_filtered$Symbol, data = df_for_aggregation, mean)

# Rename the columns for clarity
colnames(gene_symbol_means) <- c("Symbol", colnames(df_for_aggregation))

# Print the resulting DataFrame
#gene_symbol_means

# Calculate Median Expression per Gene Symbol using the aggregate function

# Remove the first column (GeneSymbol) for aggregation
df_for_aggregation <- GSE84433_expression_data_GeneSymbol_filtered[, -1]

# Calculate the median for each gene symbol

gene_symbol_medians <- aggregate(. ~ GSE84433_expression_data_GeneSymbol_filtered$Symbol, data = df_for_aggregation, median)

# Rename the columns for clarity
colnames(gene_symbol_medians) <- c("Symbol", colnames(df_for_aggregation))

GSE84433_Gene_Symbols_Median_Expression_estimate_input <- gene_symbol_medians


# dim(GSE84426_Median_Expression_estimate_input)
# [1] 25124    77

# Combine the two data.frames GSE84426 and GSE84433
GSE84426_GSE84433_Median_Expr_estimateInput <- inner_join(GSE84426_Median_Expression_estimate_input, GSE84433_Gene_Symbols_Median_Expression_estimate_input, by = "Symbol")

# dim(GSE84426_GSE84433_Median_Expr_estimateInput)
# [1] 25124   434

# Send to Excel to edit for ESTIMATE input 

write.table(GSE84433_Gene_Symbols_Median_Expression_estimate_input, "GSE84433_Gene_Symbols_Median_Expression_estimate_input.txt", sep = "\t", row.names = TRUE)

write.table(GSE84426_GSE84433_Median_Expr_estimateInput, "GSE84426_GSE84433_Median_Expr_estimateInput.txt", sep = "\t", row.names = TRUE)

# Read in the edited file for ESTIMATE input

GSE84433_Gene_Symbols_Estimate_input_Median_Expression <- read.table("GSE84433_Gene_Symbols_Estimate_input_Median_Expression.txt", header = TRUE, sep = "\t")


library(estimate)

# ESTIMATE ANALYSIS BASED ON MEDIAN EXPRESSION 
filterCommonGenes("GSE84433_Gene_Symbols_Estimate_input_Median_Expression.txt", "GSE84433_Gene_Symbols_Estimate_input_Median_Expression.gct")

filterCommonGenes("GSE84426_GSE84433_Median_Expr_estimateInput.txt", "GSE84426_GSE84433_Median_Expr_estimateInput.gct")

#[1] "Merged dataset includes 10280 genes (132 mismatched)."

estimateScore("GSE84433_Gene_Symbols_Estimate_input_Median_Expression.gct", "GSE84433_Gene_Symbols_Median_Expression_estimate_scores.gct")

estimateScore("GSE84426_GSE84433_Median_Expr_estimateInput.gct", "GSE84426_GSE84433_Median_Expr_estimate_scores.gct")

#[1] "1 gene set: StromalSignature  overlap= 140"
#[1] "2 gene set: ImmuneSignature  overlap= 141"

# ESTIMATE ANALYSIS BASED ON AVERAGE EXPRESSION

filterCommonGenes("GSE84433_Gene_Symbols_Average_Expression_estimate_input.txt", "GSE84433_Gene_Symbols_Average_Expression_estimate_input.gct")
#[1] "Merged dataset includes 10280 genes (132 mismatched)."

estimateScore("GSE84433_Gene_Symbols_Average_Expression_estimate_input.gct", "GSE84433_Gene_Symbols_Average_Expression_estimate_scores.gct")
#"1 gene set: StromalSignature  overlap= 140"
#"2 gene set: ImmuneSignature  overlap= 141"


# Transpose the estimate scores .gct file such that the samples are rows and Gene Symbols are columns
# This is necessary before merging with phenoData subset (survival info)
# Transpose using excel and rename sample column to "sample_GSM" to match phenoData subset

# Read in the estimate scores into a data.frame
GSE84426_GSE84433_Median_Expr_estimateScores <- read.table("GSE84426_GSE84433_Median_Expr_estimateScores.txt", header = TRUE, sep = "\t")

# dim(GSE84426_GSE84433_Median_Expr_estimateScores)
# [1] 433   5
# Combine phenoData for GSE84426 and GSE84433

GSE84426_GSE84433_Median_Expr_phenoData <- rbind(GSE84426_phenoData_subset, GSE84433_phenoData_subset)

# dim(GSE84426_GSE84433_Median_Expr_phenoData)
# [1] 433   7

#join phenoData with estimate scores 
GSE84433_Gene_Symbols_Median_Expression_estimateScores %>%
  inner_join(GSE84433_phenoData_subset) -> GSE84433_Median_Expression_estimateScores_phenoData

#join phenoData with estimate scores for GSE84426 and GSE84433 
GSE84426_GSE84433_Median_Expr_estimateScores %>%
  inner_join(GSE84426_GSE84433_Median_Expr_phenoData) -> GSE84426_GSE84433_Median_Expr_estimateScores_phenoData

# dim(GSE84426_GSE84433_Median_Expr_estimateScores_phenoData)
# [1] 433  11


# This doesn't work because you have to transpose the expression data first before merge with phenoData (surv info)
#GSE84433_expression_data_GeneSymbol %>%
# inner_join(GSE84433_phenoData_subset ) -> GSE84433_expression_surv_info

# GSE84433_Gene_Symbols_Estimate_input_Median_Expression data.frame has Gene Symbols as rows and samples as columns so it needs to be transposed
# such that the samples become rows and Gene Symbols become columns
# However transposing converts data.frame to array matrix so there is need to remove the V1, V2, V3 column names added

GSE84426_GSE84433_Median_Expr_estimateInput <- read.table("GSE84426_GSE84433_Median_Expr_estimateInput.txt", header = TRUE, sep = "\t")

#library(tidyverse)

# Transpose the data frame and convert it back to a data frame
GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed <- as.data.frame(t(GSE84433_Gene_Symbols_Estimate_input_Median_Expression))

GSE84426_GSE84433_Median_Expr_estimateInput_transposed <- as.data.frame(t(GSE84426_GSE84433_Median_Expr_estimateInput))

# Check dimensions of transposed data.frame 
dim(GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed)
#[1]   358 25124

# dim(GSE84426_GSE84433_Median_Expr_estimateInput_transposed)
#[1]   434 25124
# Excel has a limit of 16,384 columns and 1,048,576 rows
# Split the columns so that you get two equal files 
# While in Excel, remove row with V1, V2, V3.... Rename "Gene.Symbol" column to "sample_GSM"

GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1 <- GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed[, 1:12562]

GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part2 <- GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed[, 12563:25124]

GSE84426_GSE84433_Median_Expr_estimateInput_transposed_transposed_part1 <- GSE84426_GSE84433_Median_Expr_estimateInput_transposed[, 1:12562]

GSE84426_GSE84433_Median_Expr_estimateInput_transposed_transposed_part2 <- GSE84426_GSE84433_Median_Expr_estimateInput_transposed[, 12563:25124]

# Write the files to Excel for proper editing 
write.table(GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1, "GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1.txt", sep = "\t", row.names = TRUE)

write.table(GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part2, "GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part2.txt", sep = "\t", row.names = TRUE)

write.table(GSE84426_GSE84433_Median_Expr_estimateInput_transposed_transposed_part1, "GSE84426_GSE84433_Median_Expr_estimateInput_transposed_transposed_part1.txt", sep = "\t", row.names = TRUE)

write.table(GSE84426_GSE84433_Median_Expr_estimateInput_transposed_transposed_part2, "GSE84426_GSE84433_Median_Expr_estimateInput_transposed_transposed_part2.txt", sep = "\t", row.names = TRUE)

GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1 <- read.table("GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1.txt", header = TRUE, sep = "\t")

library(dplyr)

# First inner join between GSE84433_Median_Expression_estimateScores_phenoData and part1
result_part1 <- inner_join(GSE84433_Median_Expression_estimateScores_phenoData, GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1, by = "sample_GSM")

# Second inner join between the result from part1 and part2
GSE84433_Median_Expression_surv_info_estimateScores <- inner_join(result_part1, GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part2, by = "sample_GSM")

dim(GSE84433_Median_Expression_surv_info_estimateScores)
#[1]   357 25135

# dim(GSE84426_GSE84433_Median_Expr_surv_info_estimateScores)
# [1]   433 25135

library(survival)
library(maxstat)

# NOTE: Find out if survival time is in days or months. Presumption: It is in days. 

GSE84433_Median_Expression_surv_info_estimateScores_surv_obj <- with(GSE84433_Median_Expression_surv_info_estimateScores, Surv(as.numeric(time), as.numeric(event)))

GSE84426_GSE84433_Median_Expr_surv_info_estimateScores_surv_obj <- with(GSE84426_GSE84433_Median_Expr_surv_info_estimateScores, Surv(as.numeric(time), as.numeric(event)))


GSE84433_Median_Expr_estimateScores_maxstat_stromal_optimalcutoff <- maxstat.test(GSE84433_Median_Expression_surv_info_estimateScores_surv_obj  ~ GSE84433_Median_Expression_surv_info_estimateScores$StromalScore, data = GSE84433_Median_Expression_surv_info_estimateScores, smethod = "LogRank")

#estimated cutpoint 419.5145


GSE84433_Median_Expr_estimateScores_maxstat_IMMUNE_optimalcutoff <- maxstat.test(GSE84433_Median_Expression_surv_info_estimateScores_surv_obj  ~ GSE84433_Median_Expression_surv_info_estimateScores$ImmuneScore, data = GSE84433_Median_Expression_surv_info_estimateScores, smethod = "LogRank")
#estimated cutpoint 1900.368 

# Optimal cutoff values
cutoff_stromal <- 419.5145
cutoff_immune <- 1900.368

# Stromal and Immune cutoffs for GSE84426_GSE84433
GSE84426_GSE84433_maxstat_stromal_optimalcutoff <- maxstat.test(GSE84426_GSE84433_Median_Expr_surv_info_estimateScores_surv_obj  ~ GSE84426_GSE84433_Median_Expr_surv_info_estimateScores$StromalScore, data = GSE84426_GSE84433_Median_Expr_surv_info_estimateScores, smethod = "LogRank")
#estimated cutpoint 410.7191

GSE84426_GSE84433_maxstat_Immune_optimalcutoff <- maxstat.test(GSE84426_GSE84433_Median_Expr_surv_info_estimateScores_surv_obj  ~ GSE84426_GSE84433_Median_Expr_surv_info_estimateScores$ImmuneScore, data = GSE84426_GSE84433_Median_Expr_surv_info_estimateScores, smethod = "LogRank")
#estimated cutpoint 1877.145  

# Optimal cutoff values for GSE84426_GSE84433
cutoff_stromal <- 410.7191
cutoff_immune <- 1877.145

# Divide the immune and stromal scores into the high and low groups using optimal cutoff
GSE84426_GSE84433_surv_info_estimateScores_categorized <- GSE84426_GSE84433_Median_Expr_surv_info_estimateScores %>%
  mutate(StromalScore_group = ifelse(StromalScore > cutoff_stromal, "high", "low"),
         ImmuneScore_group = ifelse(ImmuneScore > cutoff_immune, "high", "low"))

# dim(GSE84426_GSE84433_surv_info_estimateScores_categorized)
# [1]   433 25137

# Find out size of high and low in ImmuneScore_group and StromalScore_group

table(GSE84426_GSE84433_surv_info_estimateScores_categorized$ImmuneScore_group)

# high  low 
# 100  333 

table(GSE84426_GSE84433_surv_info_estimateScores_categorized$StromalScore_group)

# high  low 
# 240  193 

# Create new columns for categorizing StromalScore and ImmuneScore
GSE84433_Median_Expression_surv_info_estimateScores$StromalScore_group <- ifelse(
  GSE84433_Median_Expression_surv_info_estimateScores$StromalScore >= cutoff_stromal,
  "high",
  "low"
)

GSE84433_Median_Expression_surv_info_estimateScores$ImmuneScore_group <- ifelse(
  GSE84433_Median_Expression_surv_info_estimateScores$ImmuneScore >= cutoff_immune,
  "high",
  "low"
)

library(dplyr)

# Rearrange the data.frame
GSE84426_GSE84433_surv_info_estimateScores_categorized <- GSE84426_GSE84433_surv_info_estimateScores_categorized %>%
  select(1:4, StromalScore_group, ImmuneScore_group, 5:ncol(.))

# Rearrange the data.frame
GSE84433_Median_Expression_surv_info_estimateScores_categorized <- GSE84433_Median_Expression_surv_info_estimateScores %>%
  select(1:4, StromalScore_group, ImmuneScore_group, 5:ncol(.))

# Convert time and event to appropriate data types
GSE84433_Median_Expression_surv_info_estimateScores_categorized$time <- as.numeric(GSE84433_Median_Expression_surv_info_estimateScores_categorized$time)
GSE84433_Median_Expression_surv_info_estimateScores_categorized$event <- as.numeric(GSE84433_Median_Expression_surv_info_estimateScores_categorized$event)

# convert ("time", "event", "age") from character to numeric using lapply function

GSE84426_GSE84433_surv_info_estimateScores_categorized[, c("time", "event", "age")] <- lapply(GSE84426_GSE84433_surv_info_estimateScores_categorized[, c("time", "event", "age")], as.numeric)


# Create a Surv object
surv_obj <- Surv(time = GSE84433_Median_Expression_surv_info_estimateScores_categorized$time, event = GSE84433_Median_Expression_surv_info_estimateScores_categorized$event)

surv_obj <- Surv(time = GSE84426_GSE84433_surv_info_estimateScores_categorized$time, event = GSE84426_GSE84433_surv_info_estimateScores_categorized$event)


#High vs Low ImmuneScore_group survival analysis 

immune_km_fit <- survfit(surv_obj ~ ImmuneScore_group, data = GSE84433_Median_Expression_surv_info_estimateScores_categorized)

GSE84426_GSE84433_immune_km_fit <- survfit(surv_obj ~ ImmuneScore_group, data = GSE84426_GSE84433_surv_info_estimateScores_categorized)

GSE84426_GSE84433_stromal_km_fit <- survfit(surv_obj ~ StromalScore_group, data = GSE84426_GSE84433_surv_info_estimateScores_categorized)

library(survminer)

survminer::ggsurvplot(immune_km_fit, data = GSE84433_Median_Expression_surv_info_estimateScores_categorized, risk.table = T, pval = TRUE)

GSE84426_GSE84433_immune_KM_SurvPlot <- ggsurvplot(GSE84426_GSE84433_immune_km_fit, data = GSE84426_GSE84433_surv_info_estimateScores_categorized, risk.table = T, pval = TRUE)

GSE84426_GSE84433_stromal_KM_SurvPlot <- ggsurvplot(GSE84426_GSE84433_stromal_km_fit, data = GSE84426_GSE84433_surv_info_estimateScores_categorized, risk.table = T, pval = TRUE)


print(GSE84426_GSE84433_immune_KM_SurvPlot)
print(GSE84426_GSE84433_stromal_KM_SurvPlot)

# ggsave(file = "GSE84426_GSE84433_immune_KM_SurvPlot", print(survp))


#High vs Low StromalScore_group survival analysis 
stromal_km_fit <- survfit(surv_obj ~ StromalScore_group, data = GSE84433_Median_Expression_surv_info_estimateScores_categorized)

survminer::ggsurvplot(stromal_km_fit, data = GSE84433_Median_Expression_surv_info_estimateScores_categorized, risk.table = T, pval = TRUE)

#count number of high vs low
GSE84433_immune_counts <- table(GSE84433_Median_Expression_surv_info_estimateScores_categorized$ImmuneScore_group)

GSE84433_stromal_counts <- table(GSE84433_Median_Expression_surv_info_estimateScores_categorized$StromalScore_group)


# Load the limma library
library(limma)

# Create a factor variable for high and low ImmuneScore

ImmuneScore_group <- factor(GSE84433_Median_Expression_surv_info_estimateScores_categorized$ImmuneScore_group, levels = c("high", "low"))

ImmuneScore_group <- factor(GSE84426_GSE84433_surv_info_estimateScores_categorized$ImmuneScore_group, levels = c("high", "low"))

# Create the design matrix for ImmuneScore_group 
immune_design_matrix <- model.matrix(~0 + ImmuneScore_group)

# Create contrasts using the immune_design_matrix
immune_contrast_matrix <- makeContrasts(high_vs_low = ImmuneScore_grouphigh - ImmuneScore_grouplow, levels = colnames(immune_design_matrix))


# Extract gene expression data from the data frame # Exclude non-gene columns

GSE84433_Median_Expression_lmFit <- GSE84433_Median_Expression_surv_info_estimateScores_categorized[, -(1:13)] 

GSE84426_GSE84433_lmFit <- GSE84426_GSE84433_surv_info_estimateScores_categorized[, -(1:13)] 


# Transpose it such that the row dimension of design matrix matches column dimension of expression data
GSE84433_Median_Expression_lmFit <- t(GSE84433_Median_Expression_lmFit)

GSE84426_GSE84433_lmFit <- t(GSE84426_GSE84433_lmFit)

# Perform differential expression analysis
immune_lmFit <- lmFit(GSE84433_Median_Expression_lmFit, immune_design_matrix)

immune_lmFit <- lmFit(GSE84426_GSE84433_lmFit, immune_design_matrix)

fit_contrast_immune <- contrasts.fit(immune_lmFit, immune_contrast_matrix)

immune_eBayes <- eBayes(fit_contrast_immune)

# Extract the differentially expressed genes

GSE84426_GSE84433_Immune_DEGs <- topTable(immune_eBayes, number = Inf)

# Get all significant Immune DEGs

GSE84426_GSE84433_Significant_Immune_DEGs_adj.P.Val_0.05 <- GSE84426_GSE84433_Immune_DEGs[GSE84426_GSE84433_Immune_DEGs$adj.P.Val < 0.05, ]

# Define the genes you want to check

genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")

# Check if the genes are present in the row names of GSE84433_Stromal_DEGS

matching_rows <- rownames(GSE84426_GSE84433_Immune_DEGs) %in% genes_to_check

# Find if genes_to_check are present
GSE84426_GSE84433_Immune_DEGs_genes_to_check <- GSE84426_GSE84433_Immune_DEGs[matching_rows, ]

# Genes are already arranged in ascending order (starting with smallest adj.P.Val)
#GSE84426_GSE84433_Immune_DEGs_genes_to_check <- GSE84426_GSE84433_Immune_DEGs_genes_to_check %>%
# arrange(adj.P.Val)


# Format the p-values in human-readable notation

GSE84426_GSE84433_Immune_DEGs_genes_to_checkFormatted <- GSE84426_GSE84433_Immune_DEGs_genes_to_check

GSE84426_GSE84433_Immune_DEGs_genes_to_checkFormatted$adj.P.Val <- format(GSE84426_GSE84433_Immune_DEGs_genes_to_checkFormatted$adj.P.Val, digits = 3, scientific=FALSE) 

# print(GSE84426_GSE84433_Immune_DEGs_genes_to_checkFormatted)


# GSE84426_GSE84433 StromalScore_group ANALYSIS 

StromalScore_group <- factor(GSE84426_GSE84433_surv_info_estimateScores_categorized$StromalScore_group, levels = c("high", "low"))

stromal_design_matrix <- model.matrix(~0 + StromalScore_group)

stromal_contrast_matrix <- makeContrasts(high_vs_low = StromalScore_grouphigh - StromalScore_grouplow, levels = colnames(stromal_design_matrix))

stromal_lmFit <- lmFit(GSE84426_GSE84433_lmFit, stromal_design_matrix)

fit_contrast_stromal <- contrasts.fit(stromal_lmFit, stromal_contrast_matrix)

stromal_eBayes <- eBayes(fit_contrast_stromal)

GSE84426_GSE84433_Stromal_DEGs <- topTable(stromal_eBayes, number = Inf)

GSE84426_GSE84433_Significant_Stromal_DEGs_adj.P.Val_0.05  <- GSE84426_GSE84433_Stromal_DEGs[GSE84426_GSE84433_Stromal_DEGs$adj.P.Val < 0.05, ]



GSE84426_GSE84433_intersection_genes <- intersect(rownames(GSE84426_GSE84433_Significant_Immune_DEGs_adj.P.Val_0.05), rownames(GSE84426_GSE84433_Significant_Stromal_DEGs_adj.P.Val_0.05))

length(GSE84426_GSE84433_intersection_genes)
# [1] 2945

# Check if the genes are present in the row names of GSE84433_Stromal_DEGS

matching_rows <- rownames(GSE84426_GSE84433_Stromal_DEGs) %in% genes_to_check

# Find if genes_to_check are present
GSE84426_GSE84433_Stromal_DEGs_genes_to_check <- GSE84426_GSE84433_Stromal_DEGs[matching_rows, ]


genes_to_check %in% GSE84426_GSE84433_intersection_genes


# Format the p-values in human-readable notation

GSE84426_GSE84433_Stromal_DEGs_genes_to_checkFormatted <- GSE84426_GSE84433_Stromal_DEGs_genes_to_check

GSE84426_GSE84433_Stromal_DEGs_genes_to_checkFormatted$adj.P.Val <- format(GSE84426_GSE84433_Stromal_DEGs_genes_to_checkFormatted$adj.P.Val, digits = 3, scientific=FALSE) 

# print(GSE84426_GSE84433_Stromal_DEGs_genes_to_checkFormatted)



length(GSE84433_intersection_genes)
# [1] 2263

# UNIVARIATE ANALYSIS 

# Remove NAs
GSE84433_univariate_analysis_input_clean_data <- GSE84433_Median_Expression_surv_info_estimateScores_categorized[!is.na(GSE84433_Median_Expression_surv_info_estimateScores_categorized$time), ]

GSE84426_GSE84433_surv_info_estimateScores_categorized_filtered <- GSE84426_GSE84433_surv_info_estimateScores_categorized[!is.na(GSE84426_GSE84433_surv_info_estimateScores_categorized$time), ]


# keep only those rows where "time" > 0
GSE84426_GSE84433_surv_info_estimateScores_categorized_431 <- GSE84426_GSE84433_surv_info_estimateScores_categorized_filtered[GSE84426_GSE84433_surv_info_estimateScores_categorized_filtered$time > 0, ]

dim(GSE84426_GSE84433_surv_info_estimateScores_categorized_431)
# [1]   431 25137

# dim(GSE84433_univariate_analysis_input_clean_data)
# [1]   355 25137


# Create the Surv() object
surv_obj <- Surv(GSE84426_GSE84433_surv_info_estimateScores_categorized_431$time, GSE84426_GSE84433_surv_info_estimateScores_categorized_431$event)

# Initialize an empty list to store results
GSE84426_GSE84433_intersection_genes_univariate_results <- list()

# Perform univariate Cox regression for each gene in GSE84426_GSE84433_intersection_genes
for (gene in GSE84426_GSE84433_intersection_genes) {
  cox_model <- coxph(surv_obj ~ get(gene), data = GSE84426_GSE84433_surv_info_estimateScores_categorized_431)
  result <- summary(cox_model)
  GSE84426_GSE84433_intersection_genes_univariate_results[[gene]] <- result
}

# Find Significant genes in GSE84426_GSE84433_intersection_genes_univariate_results


# Find intersection genes that are both in TCGA and GSE84426_GSE84433

GSE84426_GSE84433_intersection_genes %in% TCGA_Significant_intersection_genes


# Can't use all the columns, predictors too many (2263)
#GSE84433_Median_Expression_multivariate_results <- coxph(surv_obj ~ ., data = GSE84433_univariate_analysis_input_clean_data)  # Only include predictor variables


library(dplyr)

#rearranges the columns so that CCR4, FMO2, and SELP come first
GSE84433_univariate_analysis_input_clean_data <- GSE84433_univariate_analysis_input_clean_data %>%
  select(sample_GSM, StromalScore, ImmuneScore, ESTIMATEScore, 
         StromalScore_group, ImmuneScore_group, TumorPurity, age, event, time,
         Nstage, Tstage, Gender, CCR4, FMO2, SELP, everything())

# Perform multivariate Cox regression using columns 14:100

GSE84433_Median_Expression_multivariate_results <- coxph(surv_obj ~ ., data = GSE84433_univariate_analysis_input_clean_data[, 14:100])


# Perform multivariate Cox regression with SELP, CCR4, and FMO2 as predictors
multivariate_cox_model <- coxph(
  surv_obj ~ SELP + CCR4 + FMO2,
  data = GSE84433_univariate_analysis_input_clean_data
)

# Get the summary of the multivariate Cox regression
#GSE84433_Median_Expression_multivariate_results <- summary(multivariate_cox_model)
