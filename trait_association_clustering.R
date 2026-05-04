rm(list = ls())
library(data.table) # data.table_1.16.4
library(readxl) # readxl_1.4.3
library(dplyr) # dplyr_1.1.4
library(fuzzyjoin) # fuzzyjoin_0.1.6.1
library(pheatmap) # pheatmap_1.0.12
library(factoextra) # factoextra_1.0.7
library(patchwork) # patchwork_1.3.2
library(ggplot2) # ggplot2_3.5.2
library(ggplotify) # ggplotify_0.1.3 (To convert pheatmap to ggplot2)
options(width = 160)
set.seed(1)

load("Trait_Correlation_Heatmaps.RData")
dim(pet_pos_expr_nom_DAA_pheno) # 178 492

# Append age_onset related phenotypes
cutoff_age <- 65
pet_pos_expr_nom_DAA_pheno <- pet_pos_expr_nom_DAA_pheno %>%
  mutate(OnsetGroup = ifelse(Age_onset < cutoff_age, "Early", "Late"),
  OnsetGroup = factor(OnsetGroup, levels = c("Late", "Early")))
table(pet_pos_expr_nom_DAA_pheno$OnsetGroup)
# Late Early
#   75    16
pet_pos_expr_nom_DAA_pheno$DiseaseDuration <- pet_pos_expr_nom_DAA_pheno$Age_at_last - pet_pos_expr_nom_DAA_pheno$Age_onset
summary(pet_pos_expr_nom_DAA_pheno$DiseaseDuration)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#  1.000   3.500   5.000   5.824   8.000  17.000      87
  
phenotypes <- c("Zscore_AB40", "Zscore_AB42", "Zscore_AB42_AB40", "Zscore_Tau", "Zscore_pTau", "Clinical_bin", "AT_bin", "CDR_score_closest", "CDR_bin", "sumbox", "MMSE", "APOE_bin", "PC1", "PC2", "Age_onset", "OnsetGroup", "DiseaseDuration")

sumstats <- list()
for (i in 1:length(phenotypes)) {
  pheno <- phenotypes[i]
  col_to_extract <- c("UniquePhenoID", "Age_at_draw", "Sex", pheno, names(pet_pos_expr_nom_DAA_pheno)[42:ncol(pet_pos_expr_nom_DAA_pheno)])
  df <- pet_pos_expr_nom_DAA_pheno[, col_to_extract]
  df <- df[!is.na(df[,4]),]
  dim_df <- nrow(df)
  analyte <- as.data.frame(df[,5:ncol(df)])
  analyte_name <- as.data.frame(colnames(analyte))
  sumstats_df <- data.frame()
  for (j in 1:ncol(analyte)) {
    analyte_ID <- analyte_name[j,1]
    protein <- analyte[,j]
    formula <- paste0("protein ~ ", pheno, " + as.numeric(Age_at_draw) + as.factor(Sex)", sep="")
    model <- lm(formula, data = df)
    output <- as.data.frame(summary(model)[4])
    result <- cbind(as.character(analyte_ID), output[2,1], output[2,2], output[2,4])
    sumstats_df <<- rbind(sumstats_df, result)
  }
  colnames(sumstats_df) <- c("Analyte", "Estimate", "Standard_error", "Pvalue")
  cols <- c(2,3,4)
  sumstats_df[,cols] <- apply(sumstats_df[,cols], 2, function(x) as.numeric(as.character(x)))
  sumstats_df$FDR <- p.adjust(sumstats_df$Pvalue, method = "fdr", n = nrow(sumstats_df))
  sumstats_df <- sumstats_df[order(sumstats_df$Pvalue),]
  # make a variable for saving "Estimate" from each model
  variable_name <- paste0(pheno,"_Estimate")
  assign(variable_name, sumstats_df$Estimate)
  # make a dataframe for saving each model
  df_name <- paste0(pheno,"_N",dim_df)
  assign(df_name, sumstats_df)
  # assign the model dataframe to list object
  sumstats[[paste0(pheno,"_N",dim_df)]] <- sumstats_df
  print(pheno)
}
summary(sumstats)
head(Clinical_bin_N178, 2)
head(Age_onset_N91, 2)
head(z_log_odds_N178, 2)
nrow(z_log_odds_N178[z_log_odds_N178$FDR < 0.05,]) # 7
nrow(z_log_odds_N178[z_log_odds_N178$Pvalue < 0.05,]) # 72

# Get estimate for each protein against 14 phenotypes
estimate_mat <- as.data.frame(rbind(Zscore_AB40_Estimate, Zscore_AB42_Estimate, Zscore_AB42_AB40_Estimate, Zscore_Tau_Estimate, Zscore_pTau_Estimate, Clinical_bin_Estimate, 
  AT_bin_Estimate, CDR_score_closest_Estimate, CDR_bin_Estimate, sumbox_Estimate, MMSE_Estimate, Age_onset_Estimate, OnsetGroup_Estimate, DiseaseDuration_Estimate))
names(estimate_mat) <- Zscore_AB40_N110$Analyte
row.names(estimate_mat) <- c("AB40", "AB42", "AB42_AB40", "Tau", "pTau", "Clinical_status", "AT_status", "CDR", "CDR_bin", "CDRSB", "MMSE", "Age_onset", "Onset_group", "Disease_duration")
dim(estimate_mat) # 14 454 (14 phenotypes, 454 proteins)

# For protein clustering
mat_scaled_t <- t(estimate_mat)
dim(mat_scaled_t) # 454  14

# Identify optmal number of clusters
# Silhouette method
pdf("Optimal_cluster_identification_silhoutte.pdf", width = 6, height = 6)
fviz_nbclust(mat_scaled_t, FUN = hcut, method = "silhouette") # 3 clusters
dev.off()
# Gap statistic
pdf("Optimal_cluster_identification_gap_stats.pdf", width = 6, height = 6)
fviz_nbclust(mat_scaled_t, FUN = hcut, method = "gap_stat", nstart = 25, k.max = 15) # 5 clusters
dev.off()
# Total within sum of square
pdf("Optimal_cluster_identification_wss.pdf", width = 6, height = 6)
fviz_nbclust(mat_scaled_t, FUN = hcut, method = "wss") + geom_vline(xintercept = 5, linetype = 2) # 5 clusters
dev.off()

gap_stats <- fviz_nbclust(mat_scaled_t, FUN = hcut, method = "gap_stat", nstart = 25, k.max = 15)
total_wss <- fviz_nbclust(mat_scaled_t, FUN = hcut, method = "wss") + geom_vline(xintercept = 5, linetype = 2)
clustering_figure <- (gap_stats | total_wss) + plot_annotation(tag_levels = 'A')
clustering_figure
ggsave("Optimal_cluster_identification_gap_wss.png", plot = clustering_figure, width = 10, height = 5, units = "in", dpi = 600, bg = "white")
ggsave("Optimal_cluster_identification_gap_wss.pdf", plot = clustering_figure, width = 10, height = 5, units = "in", dpi = 600, bg = "white", device = cairo_pdf)

# According to Gap statistics and Total within sum of square, optimal number of clusters is 5
optimal_k_proteins <- 5

# Hierarchical clustering for proteins
d_proteins <- dist(mat_scaled_t, method = "euclidean")
hc_proteins <- hclust(d_proteins, method = "ward.D2")
# Cut tree into optimal clusters
protein_clusters <- cutree(hc_proteins, k = optimal_k_proteins)
# Hierarchical clustering for phenotypes
d_pheno <- dist(estimate_mat, method = "euclidean")
hc_pheno <- hclust(d_pheno, method = "ward.D2")

p_out <- pheatmap(estimate_mat,
         cluster_rows = hc_pheno,
         cluster_cols = hc_proteins,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         annotation_col = data.frame(ProteinCluster = factor(protein_clusters)),
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize_row = 10,
         cutree_cols = 5,
         angle_col = 45,
         main = "Phenotype × Protein clustering")
pdf("Protein_association_with_phenotype_heatmap.pdf", width = 10, height = 6)
print(p_out)
dev.off()