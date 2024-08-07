
# Merge based on submitter_id
TCGA_CESC_OS_expression_surv_info_291_488_univariate_input <- TCGA_CESC_survival_time_event_months_filtered_294_arranged %>%
  inner_join(TCGA_CESC_OS_expression_gene_read_counts_limma_484_304_Tumors_arranged, by = "submitter_id")

dim(TCGA_CESC_OS_expression_surv_info_291_488_univariate_input)
# [1] 291 488

saveRDS(TCGA_CESC_OS_expression_surv_info_291_488_univariate_input, file = "TCGA_CESC_OS_expression_surv_info_291_488_univariate_input.rds")

write.table(TCGA_CESC_OS_expression_surv_info_291_488_univariate_input, file = "TCGA_CESC_OS_expression_surv_info_291_488_univariate_input.txt", sep = "\t", row.names = TRUE)

# Convert time column to numeric 
TCGA_CESC_OS_expression_surv_info_291_488_univariate_input$time <- as.numeric(TCGA_CESC_OS_expression_surv_info_291_488_univariate_input$time)

# create surv object using TCGA_CESC_OS_expression_surv_info_291_488_univariate_input
surv_obj_TCGA_CESC_days <- Surv(time = TCGA_CESC_OS_expression_surv_info_291_488_univariate_input$time, event = TCGA_CESC_OS_expression_surv_info_291_488_univariate_input$event)

surv_obj_TCGA_CESC_days_multivariate <- Surv(time = TCGA_CESC_OS_expression_surv_info_291_65_multivariate_input$time, event = TCGA_CESC_OS_expression_surv_info_291_65_multivariate_input$event)

# Perform univariate Cox regression for each gene in TCGA_GTEx_GSE63514_GSE39001_common_genes_484

# Initialize an empty list to store results
TCGA_CESC_OS_expression_surv_info_291_488_univariate_results <- list()

for (gene in TCGA_GTEx_GSE63514_GSE39001_common_genes_484) {
  cox_model <- coxph(surv_obj_TCGA_CESC_days ~ get(gene), data = TCGA_CESC_OS_expression_surv_info_291_488_univariate_input)
  result <- summary(cox_model)
  TCGA_CESC_OS_expression_surv_info_291_488_univariate_results[[gene]] <- result
}

# Define the significance level
significance_level <- 0.05

# Extract genes significant in the univariate analysis 
TCGA_CESC_OS_univariate_significant_genes <- names(TCGA_CESC_OS_expression_surv_info_291_488_univariate_results)[sapply(TCGA_CESC_OS_expression_surv_info_291_488_univariate_results, function(result) {
  p_value <- result$coefficients["get(gene)", "Pr(>|z|)"]
  p_value < significance_level
})]

# Display the list of significant genes
print(TCGA_CESC_OS_univariate_significant_genes)

# Get the coef, Hazard Ratio/exp(coef), P-value/Pr(>|z|) and 95% CI for each gene in TCGA_CESC_OS_univariate_significant_genes
TCGA_CESC_OS_expression_surv_info_291_488_univariate_results[["TXN2"]]

# subset using TCGA_CESC_OS_univariate_significant_genes for multivariate analysis 
subsetted_data <- subset(TCGA_CESC_OS_expression_surv_info_291_488_univariate_input, select = c("submitter_id", "time", "event", TCGA_CESC_OS_univariate_significant_genes))
TCGA_CESC_OS_expression_surv_info_291_65_multivariate_input <- subsetted_data

surv_obj_TCGA_CESC_days_multivariate <- Surv(time = TCGA_CESC_OS_expression_surv_info_291_65_multivariate_input$time, event = TCGA_CESC_OS_expression_surv_info_291_65_multivariate_input$event)


# partial multivariate results for genes columns 4 to 38

TCGA_CESC_OS_expression_surv_info_291_65_multivariate_results_4_38 <- coxph(surv_obj_TCGA_CESC_days_multivariate ~ ., data = TCGA_CESC_OS_expression_surv_info_291_65_multivariate_input[, 4:38])
TCGA_CESC_OS_expression_surv_info_291_65_multivariate_results_39_68 <- coxph(surv_obj_TCGA_CESC_days_multivariate ~ ., data = TCGA_CESC_OS_expression_surv_info_291_65_multivariate_input[, 39:68])


TCGA_CESC_OS_expression_surv_info_291_488_multivariate_results
Call: coxph(formula = surv_obj_TCGA_CESC_days ~ ., data = TCGA_CESC_OS_expression_surv_info_291_488_univariate_input[, 455:488])



