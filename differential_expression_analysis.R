rm(list = ls())
options(stringsAsFactors = FALSE)

library(readr)
library(readxl)
library(data.table)
library(dplyr)
library(broom)
library(writexl)
library(patchwork)
library(ggplot2)
library(ggrepel)

data <- read_excel("cleaned_data_558samples_6873proteins.xlsx")

names(data)[1:20]
table(data$final_decision)
mean(data$Age_at_draw, na.rm = TRUE)
sd(data$Age_at_draw, na.rm = TRUE)

protein_cols <- names(data)[
  grepl("^X", names(data)) &
    sapply(data, is.numeric)
]

run_lm_by_protein <- function(data_input, protein_cols, predictor, term_name) {
  
  results <- lapply(protein_cols, function(protein) {
    
    formula <- as.formula(paste0(
      protein,
      " ~ as.numeric(", predictor, ") + as.numeric(Age_at_draw) + Sex + as.numeric(PC1) + as.numeric(PC2) + as.factor(Tracer)"
    ))
    
    model <- lm(formula, data = data_input)
    
    broom::tidy(model) %>%
      dplyr::filter(term == term_name) %>%
      dplyr::mutate(Protein = protein)
  })
  
  final_results <- dplyr::bind_rows(results) %>%
    dplyr::select(Protein, estimate, std.error, statistic, p.value) %>%
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH"))
  
  return(final_results)
}

###########################
#  Dichotomized analysis  #
###########################

data <- data %>%
  dplyr::mutate(
    Biomarker_bin = ifelse(
      final_decision == "negative", 1,
      ifelse(final_decision == "positive", 2, NA)
    )
  )

table(data$Biomarker_bin, useNA = "ifany")
# 1   2 
# 380 178 

Dichotomized_results <- run_lm_by_protein(
  data_input = data,
  protein_cols = protein_cols,
  predictor = "Biomarker_bin",
  term_name = "as.numeric(Biomarker_bin)"
)

Dichotomized_results <- Dichotomized_results %>%
  dplyr::select(-statistic)

cat("Nominal p < 0.05 (dichotomized):", sum(Dichotomized_results$p.value < 0.05), "\n") # 895
cat("FDR < 0.05 (dichotomized):", sum(Dichotomized_results$FDR < 0.05), "\n") # 53

write_xlsx(
  Dichotomized_results,
  "Dichotomized_results_Tracer_895nominal_53FDR.xlsx"
)

###########################
#  Continuous: All samples #
###########################

Continuous_ALL_results <- run_lm_by_protein(
  data_input = data,
  protein_cols = protein_cols,
  predictor = "combined_amyloid_zscore",
  term_name = "as.numeric(combined_amyloid_zscore)"
)

cat("Nominal p < 0.05 (continuous all):", sum(Continuous_ALL_results$p.value < 0.05), "\n") #470
cat("FDR < 0.05 (continuous all):", sum(Continuous_ALL_results$FDR < 0.05), "\n") # 14

write_xlsx(
  Continuous_ALL_results,
  "Continuous_ALL_results_470nominal_14FDR.xlsx"
)

###########################
#  Continuous: Positive only #
###########################

data_pos <- data %>%
  dplyr::filter(final_decision == "positive")

dim(data_pos)

protein_cols_pos <- names(data_pos)[
  grepl("^X", names(data_pos)) &
    sapply(data_pos, is.numeric)
]

Continuous_POS_results <- run_lm_by_protein(
  data_input = data_pos,
  protein_cols = protein_cols_pos,
  predictor = "combined_amyloid_zscore",
  term_name = "as.numeric(combined_amyloid_zscore)"
)

cat("Nominal p < 0.05 (positive only):", sum(Continuous_POS_results$p.value < 0.05), "\n")
cat("FDR < 0.05 (positive only):", sum(Continuous_POS_results$FDR < 0.05), "\n")

write_xlsx(
  Continuous_POS_results,
  "Continuous_POS_results_454nominal_0FDR.xlsx"
)