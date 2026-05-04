rm(list = ls())
library(data.table) # data.table_1.17.8
library(readxl) # readxl_1.4.5
library(dplyr) # dplyr_1.1.4
library(fuzzyjoin) # fuzzyjoin_0.1.6.1
library(pheatmap) # pheatmap_1.0.13
library(ggplot2) # ggplot2_4.0.1
library(tidyr) # tidyr_1.3.1
library(broom) # broom_1.0.10
library(ggpubr) # ggpubr_0.6.2
library(forcats) # forcats_1.0.1
library(patchwork) # patchwork_1.3.2
options(width = 160)
set.seed(1)

# read complete expression matrix
expr_df <- data.table(read_excel("cleaned_data_558samples_6873proteins.xlsx", sheet = "Sheet1"))
dim(expr_df) #  558 6885
expr_df[1:4,1:15]

# read phenotype file
clinical_pheno <- data.table(read_excel("Phenotype_CSF_Plasma_20260204.xlsx", sheet = "Sheet1"))
dim(clinical_pheno) # 9775   30

# create expr_pheno data
expr_clin_df <- inner_join(clinical_pheno, expr_df, by=c("UniquePhenoID","DrawDate", "Age_at_draw", "Sex", "PET_Date"))
dim(expr_clin_df) # 558 6910
expr_clin_df[1:2,1:38]
table(expr_clin_df$Simplified_Status_at_draw)
#  AD ADAD   CO  CVD  DLB  FTD   OT VCID
#  82    5  441    1    4    2   21    2
table(expr_clin_df$CDR_score_closest)
#   0 0.5   1
# 447  81  30
table(expr_clin_df$CDR_score_closest, expr_clin_df$final_decision)
#      negative positive
#  0        350       97
#  0.5       25       56
#  1          5       25

cols_to_extract <- c("UniquePhenoID", "DrawDate", "APOE", "Age_at_draw", "Age_onset", "Age_at_last", "Sex", 
  "Simplified_Final_Status", "Simplified_Status_at_draw", "CDR_score_closest", "Last_CDR_score", "CSF_AT_class", "PET_Date", "tracer",
  "combined_amyloid_zscore", "final_decision", "PC1", "PC2", colnames(expr_clin_df)[38:ncol(expr_clin_df)]) 
length(cols_to_extract) # 6891
expr_clin_df <- as.data.frame(expr_clin_df)
expr_clin_df <- expr_clin_df[, cols_to_extract]
dim(expr_clin_df) # 558 6891

# get CDR-SB and MMSE scores from Knight_ADRC pheno file
MAP_pheno <- data.table(read_excel("b4_cdr.xlsx", sheet = "b4_cdr"))
dim(MAP_pheno) # 20314    29
MAP_pheno$UniquePhenoID <- paste0("MAP_", MAP_pheno$ID, sep="")
names(MAP_pheno)[2] <- "DrawDate"
MAP_pheno <- MAP_pheno[, c("UniquePhenoID", "DrawDate", "MMSE", "sumbox")]
table(expr_clin_df$UniquePhenoID %in% MAP_pheno$UniquePhenoID)
# FALSE  TRUE
#    19   539
temp <- inner_join(expr_clin_df, MAP_pheno, by=c("UniquePhenoID", "DrawDate"))
dim(temp) #  498 6893
# Since there is not a single sample with same UniquePhenoID and DrawDate, we will use 6 months-window to match.
matched_all <- fuzzy_inner_join(
  expr_clin_df %>% select(UniquePhenoID, DrawDate),
  MAP_pheno %>% select(UniquePhenoID, DrawDate),
  by = c("UniquePhenoID" = "UniquePhenoID",
         "DrawDate" = "DrawDate"),
  match_fun = list(`==`, function(x, y) abs(as.numeric(difftime(x, y, units = "days"))) <= 183)
)
matched_all <- data.table(matched_all)
dim(matched_all) # 538   4
names(matched_all)[3:4]<- c("UniquePhenoID", "DrawDate")
matched_all <- matched_all[, c(3:4)]
# make MAP_pheno cross-sectional and keep only samples within 6-month window to amyloid-PET data
MAP_pheno_6m <- inner_join(MAP_pheno, matched_all, by=c("UniquePhenoID", "DrawDate"))
MAP_pheno_6m$DrawDate <- NULL # Because DrawDate is not being used for merge and is not same as expression DrawDate
dim(MAP_pheno_6m) # 538   3

# merge CDR-SB and MMSE with actual pheno file
expr_clin_df <- left_join(expr_clin_df, MAP_pheno_6m, by="UniquePhenoID")
dim(expr_clin_df) # 558 6893
expr_clin_df <- expr_clin_df %>% select(UniquePhenoID, DrawDate, MMSE, sumbox, everything())
dim(expr_clin_df) # 558 6893
colnames(expr_clin_df)[4] <- "CDRSB"
colnames(expr_clin_df)[16:18] <- c("PET_Tracer", "PET_Zscore", "PET_Status")
expr_clin_df[1:2,1:25]


## Get Sumstats for PET asociations
# read DAA results for Continuous amyloidPET trait
pet_pos_DAA <- data.table(read_excel("Continuous_POS_results_454nominal_0FDR.xlsx", sheet = "Sheet1"))
dim(pet_pos_DAA) # 6873    5
colnames(pet_pos_DAA) <- c("Protein", "Estimate", "SE", "P", "FDR")
nrow(pet_pos_DAA[pet_pos_DAA$FDR < 0.05,]) # 0 FDR proteins
pet_pos_DAA_nom <- pet_pos_DAA[pet_pos_DAA$P < 0.05,]
dim(pet_pos_DAA_nom) # 454   5

# read soma7k annotation file
annot <- read.table("Plasma_SOMAscan7k_analyte_information.tsv", sep="\t", header=T, stringsAsFactors=F, quote="")
dim(annot) # 7291   12
annot <- annot[,c("Analytes", "EntrezGeneSymbol")]
names(annot) <- c("Protein", "Symbol")
pet_pos_DAA_nom <- inner_join(pet_pos_DAA_nom, annot, by="Protein")
dim(pet_pos_DAA_nom) # 454   6
nrow(pet_pos_DAA_nom[pet_pos_DAA_nom$FDR < 0.05,]) # 0 FDR proteins

# Bio-Hermes AB+
external_pos_DAA <- data.table(read_excel("Bio-Hermes_centiloid_Sumstats.xlsx", sheet = "continuous_POS"))
external_pos_DAA <- external_pos_DAA %>% select(Protein, estimate, std.error, p.value, FDR)
colnames(external_pos_DAA) <- c("Protein", "Estimate", "SE", "P", "FDR")
nrow(external_pos_DAA[external_pos_DAA$FDR < 0.05,]) # 10 FDR proteins
nrow(external_pos_DAA[external_pos_DAA$P < 0.05,]) # 886 Nominal proteins
external_pos_DAA <- external_pos_DAA[order(external_pos_DAA$P),]
dim(external_pos_DAA) # 6665    5
external_pos_DAA_nom <- external_pos_DAA[external_pos_DAA$P < 0.05,]
dim(external_pos_DAA_nom) # 886   5

# Overlap of nominal proteins in discovery
overlapping <- inner_join(pet_pos_DAA_nom, external_pos_DAA, by="Protein")
dim(overlapping) # 436  10

concordant <- rbind(overlapping[overlapping$Estimate.x > 0 & overlapping$Estimate.y > 0,], 
  overlapping[overlapping$Estimate.x < 0 & overlapping$Estimate.y < 0,])
nrow(concordant) # 290
concordant_nom <- as.data.frame(concordant[concordant$P.y < 0.05,])
dim(concordant_nom) # 54 10
colnames(concordant_nom) <- c("Protein", "B_KADRC", "SE_KADRC", "P_KADRC", "FDR_KADRC", "Symbol", "B_BH", "SE_BH", "P_BH", "FDR_BH")


#################################
## Proteomic Score Calculation ##
#################################

# STEP 1: Prepare protein list and weights
# Extract protein list
protein_list <- concordant_nom$Protein
# Extract weights (Knight-ADRC betas)
weights <- concordant_nom$B_KADRC
# Name weights by protein
names(weights) <- protein_list
# Confirm proteins exist in expression dataframe
missing_proteins <- setdiff(protein_list, colnames(expr_clin_df))
print(missing_proteins)

# STEP 2: Z-score standardization of proteins (Already done)
# Create copy
df <- expr_clin_df

# STEP 3: Construct weighted proteomic score
# Convert to matrix
protein_matrix <- as.matrix(expr_clin_df[, protein_list])
# Ensure order matches weights
protein_matrix <- protein_matrix[, names(weights)]
# Weight vector
w <- weights
# Create matrix of weights repeated per individual
weight_matrix <- matrix(
  rep(w, each = nrow(protein_matrix)),
  nrow = nrow(protein_matrix)
)
# Set weights to 0 where protein is missing
weight_matrix[is.na(protein_matrix)] <- 0
# Set missing protein values to 0 so they don't contribute
protein_matrix[is.na(protein_matrix)] <- 0
# Compute weighted sum
weighted_sum <- rowSums(protein_matrix * weight_matrix)
# Compute sum of absolute weights used
weight_sum <- rowSums(abs(weight_matrix))
# Compute normalized score
df$proteomic_score <- weighted_sum / weight_sum
# Z-score final score
df$proteomic_score_z <- scale(df$proteomic_score)
# Check missing
summary(df$proteomic_score_z)
#       V1
# Min.   :-3.02039
# 1st Qu.:-0.66093
# Median : 0.02091
# Mean   : 0.00000
# 3rd Qu.: 0.64695
# Max.   : 3.41676

# STEP 4: Verify score distribution
ggplot(df, aes(x = proteomic_score_z)) +
  geom_histogram(bins = 40) +
  theme_classic() +
  labs(title = "Distribution of amyloid-associated proteomic score")

# STEP 5: Association with cognitive outcomes
model_mmse <- lm(
  MMSE ~ proteomic_score_z + Age_at_draw + Sex + APOE,
  data = df
)
summary(model_mmse)
tidy(model_mmse)
# estimate = -0.0568, p = 0.57641
model_cdrsb <- lm(
  CDRSB ~ proteomic_score_z + Age_at_draw + Sex + APOE,
  data = df
)
summary(model_cdrsb)
# estimate = 0.056196, p = 0.364673
model_pet <- lm(
  PET_Zscore ~ proteomic_score_z + Age_at_draw + Sex + APOE,
  data = df
)
summary(model_pet)
# estimate = 0.037589, p = 0.3886

# STEP 6: Create clinical group variable
df$clinical_group <- df$Simplified_Status_at_draw

# STEP 7: FIGURE PANEL A
# Boxplot across clinical groups
df$clinical_group <- factor(df$clinical_group, levels=c("CO", "AD", "ADAD", "CVD", "DLB", "FTD", "OT", "VCID"))
comparisons <- list(c("CO", "AD"))
boxplot_clinical <- ggplot(df[df$clinical_group %in% c("AD", "CO") & df$PET_Status == "positive",], aes(x = clinical_group, y = proteomic_score_z, fill = clinical_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = clinical_group), width = 0.15, alpha = 0.2) +
  scale_fill_manual(values = c("AD" = "#994F00", "CO" = "#1A85FF")) +
  scale_color_manual(values = c("AD" = "#994F00", "CO" = "#1A85FF")) +
 # stat_compare_means(method = "wilcox.test", label.y = max(df[df$clinical_group %in% c("AD", "CO") & df$PET_Status == "positive",]$proteomic_score_z) * 1.1) +
  theme_classic(base_size = 14) + theme(legend.position = "none") + stat_compare_means(comparisons=comparisons) +
  labs(
    x = "Clinical group",
    y = "Proteomic score (z)",
  #  title = "Proteomic score across clinical groups"
  )
boxplot_clinical
# Optional statistical test:
compare_means(
  proteomic_score_z ~ clinical_group,
  data = df[df$clinical_group %in% c("AD", "CO") & df$PET_Status == "positive",]
)
# A tibble: 1 × 8
#  .y.               group1 group2      p p.adj p.format p.signif method
#  <chr>             <chr>  <chr>   <dbl> <dbl> <chr>    <chr>    <chr>
#1 proteomic_score_z CO     AD     0.0218 0.022 0.022    *        Wilcoxon
df$CDR_Global <- as.factor(df$CDR_score_closest)
# comparisons <- list(c("0", "0.5"), c("0.5", "1"), c("0", "1"))
comparisons <- list(c("0", "1"))
boxplot_CDR <- ggplot(df[df$PET_Status == "positive",], aes(x = CDR_Global, y = proteomic_score_z, fill = CDR_Global)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = CDR_Global), width = 0.15, alpha = 0.2) +
  scale_fill_manual(values = c("1" = "#994F00", "0" = "#1A85FF", "0.5" = "#E1BE6A")) +
  scale_color_manual(values = c("1" = "#994F00", "0" = "#1A85FF", "0.5" = "#E1BE6A")) +
  theme_classic(base_size = 14) + theme(legend.position = "none") + stat_compare_means(comparisons=comparisons) +
  labs(
    x = "Global CDR",
    y = "Proteomic score (z)",
  #  title = "Proteomic score across CDR"
  )
boxplot_CDR
# Optional statistical test:
compare_means(
  proteomic_score_z ~ CDR_Global,
  data = df[df$PET_Status == "positive",]
)
# A tibble: 3 × 8
#  .y.               group1 group2      p p.adj p.format p.signif method
#  <chr>             <chr>  <chr>   <dbl> <dbl> <chr>    <chr>    <chr>
#1 proteomic_score_z 0      0.5    0.167  0.33  0.167    ns       Wilcoxon
#2 proteomic_score_z 0      1      0.0220 0.066 0.022    *        Wilcoxon
#3 proteomic_score_z 0.5    1      0.192  0.33  0.192    ns       Wilcoxon

df$CSF_AT_class <- as.factor(df$CSF_AT_class)
comparisons <- list(c("A-T-", "A+T+"), c("A-T+", "A+T+"), c("A+T-", "A+T+"))
boxplot_AT <- ggplot(df[df$PET_Status == "positive" & df$CSF_AT_class %in% c("A-T-","A-T+","A+T-","A+T+"),], aes(x = CSF_AT_class, y = proteomic_score_z, , fill = CSF_AT_class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = CSF_AT_class), width = 0.15, alpha = 0.2) +
  scale_fill_manual(values = c("A+T+" = "#994F00", "A-T-" = "#1A85FF", "A+T-" = "#FFC20A", "A-T+" = "#E1BE6A")) +
  scale_color_manual(values = c("A+T+" = "#994F00", "A-T-" = "#1A85FF", "A+T-" = "#FFC20A", "A-T+" = "#E1BE6A")) +
  theme_classic(base_size = 14) + theme(legend.position = "none") + stat_compare_means(comparisons=comparisons) +
  labs(
    x = "CSF AT status",
    y = "Proteomic score (z)",
  #  title = "Proteomic score across CSF AT"
  )
boxplot_AT
# Optional statistical test:
compare_means(
  proteomic_score_z ~ CSF_AT_class,
  data = df[df$PET_Status == "positive" & df$CSF_AT_class %in% c("A-T-","A-T+","A+T-","A+T+"),]
)

# STEP 8: FIGURE PANEL B
# Scatter plot with MMSE
Regression_MMSE <- ggplot(df[df$PET_Status == "positive",], aes(x = proteomic_score_z, y = MMSE)) +
  geom_point(color = "#E1BE6A", alpha = 0.2) +
  geom_smooth(method = "lm", color = "#994F00", fill = "#E1BE6A", alpha = 0.4) +
 # stat_regline_equation(label.x = -2, label.y = 30) + # Regression coefficient
  stat_cor(method = "pearson", label.x = -3, label.y = 10) +
  theme_classic(base_size = 14) +
  labs(
    x = "Proteomic score (z)",
    y = "MMSE",
  #  title = "Association between proteomic score and cognition"
  )
Regression_MMSE
# Add correlation:
cor.test(df[df$PET_Status == "positive",]$proteomic_score_z, df[df$PET_Status == "positive",]$MMSE, use="complete.obs") # cor = -0.1323316, p-value = 0.08086
Regression_CDR <- ggplot(df[df$PET_Status == "positive",], aes(x = proteomic_score_z, y = CDR_score_closest)) +
  geom_point(color = "#E1BE6A", alpha = 0.2) +
  geom_smooth(method = "lm", color = "#994F00", fill = "#E1BE6A", alpha = 0.4) +
  # stat_regline_equation(label.x = -2, label.y = 30) + # Regression coefficient
  stat_cor(method = "pearson", label.x = -3, label.y = 0.75) +
  theme_classic(base_size = 14) +
  labs(
    x = "Proteomic score (z)",
    y = "CDR",
#    title = "Association between proteomic score and cognition"
  )
Regression_CDR
# Add correlation:
cor.test(df[df$PET_Status == "positive",]$proteomic_score_z, df[df$PET_Status == "positive",]$CDR_score_closest, use="complete.obs") # cor = 0.1567909, p-value = 0.03661
Regression_CDRSB <- ggplot(df[df$PET_Status == "positive",], aes(x = proteomic_score_z, y = CDRSB)) +
  geom_point(color = "#E1BE6A", alpha = 0.2) +
  geom_smooth(method = "lm", color = "#994F00", fill = "#E1BE6A", alpha = 0.4) +
  # stat_regline_equation(label.x = -2, label.y = 30) + # Regression coefficient
  stat_cor(method = "pearson", label.x = -3, label.y = 7.5) +
  theme_classic(base_size = 14) +
  labs(
    x = "Proteomic score (z)",
    y = "CDR-SB",
#    title = "Association between proteomic score and cognition"
  )
Regression_CDRSB
# Add correlation:
cor.test(df[df$PET_Status == "positive",]$proteomic_score_z, df[df$PET_Status == "positive",]$CDRSB, use="complete.obs") # cor = 0.1850926, p-value = 0.0142

# STEP 9: FIGURE PANEL C
# Forest plot of regression effects
# Test multiple outcomes:
# df$PET_bin <- ifelse(df$PET_Status == "positive", 1, 0) # This will return no result when we will focus on AB+ individuals only as it'll be just 1 group
df$Clinical_bin <- ifelse(df$clinical_group == "AD", 1, 0)
outcomes <- c(
  "MMSE",
  "CDRSB",
  "CDR_score_closest",
  "PET_Zscore",
  "Clinical_bin"
)
results <- lapply(outcomes, function(outcome){
  model <- lm(
    as.formula(
      paste(outcome, "~ proteomic_score_z + Age_at_draw + Sex + APOE + PC1 + PC2")
    ),
    #data = df
    data = df[df$PET_Status == "positive",]
  )
  tidy(model) %>%
    filter(term == "proteomic_score_z") %>%
    mutate(outcome = outcome)
}) %>%
  bind_rows()
results$outcome <- ifelse(results$outcome == "CDRSB", "CDR-SB", results$outcome)
results$outcome <- ifelse(results$outcome == "CDR_score_closest", "CDR", results$outcome)
results$outcome <- ifelse(results$outcome == "Clinical_bin", "Clinical Status", results$outcome)
results$outcome <- ifelse(results$outcome == "PET_Zscore", "Amyloid PET", results$outcome)

p3 <- ggplot(results,
             aes(x = estimate,
                 y = fct_reorder(outcome, estimate), color = outcome)) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(
      xmin = estimate - 1.96 * std.error,
      xmax = estimate + 1.96 * std.error
    ),
    height = 0.2
  ) +
  scale_color_manual(values = c(
    "MMSE"            = "#1b9e77",
    "CDR-SB"          = "#d95f02",
    "CDR"             = "#7570b3",
    "Amyloid PET"     = "#e7298a",
    "Clinical Status" = "#66a61e"
    )) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 14) + theme(legend.position = "none") +
  labs(
    x = "Effect size (beta)",
    y = "Outcome",
 #   title = "Association of proteomic score with clinical measures"
  )
p3

# STEP 10: Combine panels into one publication figure
library(patchwork)
final_figure <- (boxplot_clinical | boxplot_CDR | boxplot_AT) / (Regression_MMSE | Regression_CDRSB | p3) + 
plot_annotation(tag_levels = 'A') +
plot_layout(heights = c(2, 1.75))
final_figure

# STEP 11: Export figures (high resolution)
ggsave("proteomic_score_associations_with_traits.pdf", plot = final_figure)
ggsave("proteomic_score_associations_with_traits.png", final_figure, width = 12, height = 9, dpi = 600)

# STEP 12: Export regression results for manuscript table
write.csv(
  results,
  "proteomic_score_association_results.csv",
  row.names = FALSE
)

# Add APOE-stratified plot:
df$APOE4 <- ifelse(df$APOE %in% c(24, 34, 44), "Positive", "Negative")
APOE_MMSE <- ggplot(df, aes(x = proteomic_score_z, y = MMSE, color = factor(APOE4))) +
  geom_point() +
  geom_smooth(method = "lm") + theme_classic(base_size = 14)
APOE_MMSE <- ggplot(df[df$PET_Status == "positive",],
       aes(x = proteomic_score_z, y = MMSE, color = factor(APOE4))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(aes(color = factor(APOE4)),
           method = "pearson",
           label.x = -3,
           label.y = c(12, 10),
           show.legend = FALSE) +
  theme_classic(base_size = 14) + labs(x = "Proteomic score (z)", y = "MMSE", color = "APOE4 Status") + theme(legend.position = "none")

APOE_CDRSB <- ggplot(df, aes(x = proteomic_score_z, y = CDRSB, color = factor(APOE4))) +
  geom_point() +
  geom_smooth(method = "lm") + theme_classic(base_size = 14)
APOE_CDRSB <- ggplot(df[df$PET_Status == "positive",],
       aes(x = proteomic_score_z, y = CDRSB, color = factor(APOE4))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(aes(color = factor(APOE4)),
           method = "pearson",
           label.x = -3,
           label.y = c(10, 9),
           show.legend = FALSE) +
  theme_classic(base_size = 14) + labs(x = "Proteomic score (z)", y = "CDR-SB", color = "APOE4 Status")

APOE_figure <- (APOE_MMSE | APOE_CDRSB) + 
plot_annotation(tag_levels = 'A')
APOE_figure

# STEP 11: Export figures (high resolution)
ggsave("proteomic_score_associations_with_traits_APOE.pdf", plot = APOE_figure)  
ggsave("proteomic_score_associations_with_traits_APOE.png", APOE_figure, width = 10, height = 4, dpi = 600)

save(boxplot_clinical, boxplot_CDR, boxplot_AT, Regression_MMSE, Regression_CDRSB, p3, APOE_MMSE, APOE_CDRSB, file="proteomic_score_associations_with_traits.RData")