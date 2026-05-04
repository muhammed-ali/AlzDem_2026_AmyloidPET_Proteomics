rm(list = ls())
library(data.table) # data.table_1.16.4
library(readxl) # readxl_1.4.3
library(dplyr) # dplyr_1.1.4
library(ggplot2) # ggplot2_3.5.2
library(tidyr) # tidyr_1.3.1
library(clusterProfiler) # clusterProfiler_4.12.6
library(enrichplot) # enrichplot_1.24.4
library(org.Hs.eg.db) # org.Hs.eg.db_3.19.1
library(ReactomePA) # v1.48.0
library(DOSE) # v3.30.5
options(width = 160)
set.seed(1)

# read DAA results for amyloidPET+ individuals only
pet_pos_DAA <- data.table(read_excel("Continuous_POS_results_454nominal_0FDR.xlsx", sheet = "Sheet1"))
dim(pet_pos_DAA) # 6873    5
pet_pos_DAA_nom <- pet_pos_DAA[pet_pos_DAA$p.value < 0.05,]
dim(pet_pos_DAA_nom) # 454   5
# read soma7k annotation file
annot <- read.table("Plasma_SOMAscan7k_analyte_information.tsv", sep="\t", header=T, stringsAsFactors=F, quote="")
dim(annot) # 7291   12
annot <- annot[,c("Analytes", "UniProt", "EntrezGeneID", "EntrezGeneSymbol")]
colnames(annot)[1] <- "Protein"

pet_pos_DAA_annot <- inner_join(pet_pos_DAA, annot, by="Protein")
dim(pet_pos_DAA_annot) # 6873    8
pet_pos_DAA_nom_annot <- inner_join(pet_pos_DAA_nom, annot, by="Protein")
dim(pet_pos_DAA_nom_annot) # 454   8

nom_annot_split <- pet_pos_DAA_nom_annot %>%
  separate_rows(UniProt, EntrezGeneID, EntrezGeneSymbol, sep = "\\|")
nom_annot_split <- nom_annot_split[order(nom_annot_split$p.value),]
length(unique(nom_annot_split$EntrezGeneSymbol)) # 451
length(unique(nom_annot_split$UniProt)) # 451
length(unique(nom_annot_split$EntrezGeneID)) # 451
head(nom_annot_split[duplicated(nom_annot_split$EntrezGeneSymbol),]) # view duplicated rows
nom_annot_split <- nom_annot_split[!duplicated(nom_annot_split$EntrezGeneSymbol),]
dim(nom_annot_split) # 451   8

comp_annot_split <- pet_pos_DAA_annot %>%
  separate_rows(UniProt, EntrezGeneID, EntrezGeneSymbol, sep = "\\|")
comp_annot_split <- comp_annot_split[order(comp_annot_split$p.value),]
length(unique(comp_annot_split$EntrezGeneSymbol)) # 6049
length(unique(comp_annot_split$UniProt)) # 6066
length(unique(comp_annot_split$EntrezGeneID)) # 6056
head(comp_annot_split[duplicated(comp_annot_split$EntrezGeneSymbol),]) # view duplicated rows
comp_annot_split <- comp_annot_split[!duplicated(comp_annot_split$EntrezGeneSymbol),]
dim(comp_annot_split) # 6049    8

# load clustered object
load("Trait_Association_Heatmap.RData")
lengths(cluster_list)
#  1   2   3   4   5
# 84  51  59  80 180
head(cluster_list[[1]])
# Convert the list into a long lookup table
cluster_lookup <- stack(cluster_list) %>%
  dplyr::rename(Protein = values, Cluster = ind) %>%
  mutate(Cluster = paste0("C", Cluster))
# Join back to your dataframe
nom_annot_split <- nom_annot_split %>%
  inner_join(cluster_lookup, by = "Protein")
dim(nom_annot_split) # 451   9
table(nom_annot_split$Cluster)
# C1  C2  C3  C4  C5
# 83  51  57  78 182

# convert EntrezGeneID column to numeric (some may be characters)
nom_annot_split$EntrezGeneID <- as.character(nom_annot_split$EntrezGeneID)
comp_annot_split$EntrezGeneID <- as.character(comp_annot_split$EntrezGeneID)

save(nom_annot_split, comp_annot_split, file="nom_comp_annot_split.RData")


##############################
## Gene Ontology (EnrichGO) ##
##############################

# Run enrichGO for each cluster
enrich_results <- lapply(unique(nom_annot_split$Cluster), function(cl){
  genes <- na.omit(nom_annot_split$EntrezGeneID[nom_annot_split$Cluster == cl])
  eg <- enrichGO(gene = genes,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pAdjustMethod = "fdr", 
                 readable = TRUE)
  if(!is.null(eg)){
    as.data.frame(eg) %>% mutate(Cluster = cl)
  }
})

# Combine all clusters
enrich_df <- bind_rows(enrich_results)
table(enrich_df$Cluster, enrich_df$ONTOLOGY)
#      BP  CC  MF
#  C1   2   0   0
#  C4   0   6   0
#  C5 257   4  12

# Add enrichment score (fold enrichment or gene ratio)
enrich_df <- enrich_df %>%
  mutate(negLogFDR = -log10(p.adjust))

enrich_top10 <- enrich_df %>%
  group_by(Cluster, ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup()
dim(enrich_top10) # 21 15
table(enrich_top10$Cluster, enrich_top10$ONTOLOGY)
#     BP CC MF
#  C1  2  0  0
#  C4  0  5  0
#  C5  5  4  5

# Ensure ONTOLOGY is factor with fixed order
enrich_top10$ONTOLOGY <- paste0("GO:", enrich_top10$ONTOLOGY, sep="")
enrich_top10$ONTOLOGY <- factor(enrich_top10$ONTOLOGY, levels = c("GO:BP", "GO:CC", "GO:MF"))

# Plot
ggplot(enrich_top10, aes(x = FoldEnrichment,
                      y = reorder(Description, FoldEnrichment),
                      size = negLogFDR,
                      color = Cluster)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_size_continuous(name = expression(-log[10]~FDR)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  xlab("Fold enrichment") + ylab("Pathways")


##############################
##### KEGG (enrichKEGG) ######
##############################

# Run enrichKEGG for each cluster
enrich_results_KEGG <- lapply(unique(nom_annot_split$Cluster), function(cl){
  genes <- na.omit(nom_annot_split$EntrezGeneID[nom_annot_split$Cluster == cl])
  eg <- enrichKEGG(gene = genes,
                 organism = "hsa",
                 pAdjustMethod = "fdr", 
                 use_internal_data = FALSE)
  if(!is.null(eg)){
    as.data.frame(eg) %>% mutate(Cluster = cl)
  }
})

# Combine all clusters
enrich_KEGG <- bind_rows(enrich_results_KEGG)
table(enrich_KEGG$Cluster)
# C5
# 21

# Add enrichment score (fold enrichment or gene ratio)
enrich_KEGG <- enrich_KEGG %>%
  mutate(negLogFDR = -log10(p.adjust))

enrichKEGG_top10 <- enrich_KEGG %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup()
dim(enrichKEGG_top10) # 5 15
table(enrichKEGG_top10$Cluster)
# C5
#  5


# Ensure ONTOLOGY is factor with fixed order
enrichKEGG_top10$ONTOLOGY <- "KEGG"
# Plot
ggplot(enrichKEGG_top10, aes(x = FoldEnrichment,
                      y = reorder(Description, FoldEnrichment),
                      size = negLogFDR,
                      color = Cluster)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_size_continuous(name = expression(-log[10]~FDR)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  xlab("Fold enrichment") + ylab("Pathways")


##############################
## Reactome (enrichPathway) ##
##############################

# Run enrichPathway for each cluster
enrich_results_RCTM <- lapply(unique(nom_annot_split$Cluster), function(cl){
  genes <- na.omit(nom_annot_split$EntrezGeneID[nom_annot_split$Cluster == cl])
  eg <- enrichPathway(gene = genes,
                 organism = "human",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr", 
                 readable = TRUE)
  if(!is.null(eg)){
    as.data.frame(eg) %>% mutate(Cluster = cl)
  }
})
# Combine all clusters
enrich_RCTM <- bind_rows(enrich_results_RCTM)
table(enrich_RCTM$Cluster)
# C2 C5
#  9  5

# Add enrichment score (fold enrichment) and significance
enrich_RCTM <- enrich_RCTM %>%
  mutate(FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio),
    negLogFDR = -log10(p.adjust))

enrichRCTM_top10 <- enrich_RCTM %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup()
dim(enrichRCTM_top10) # 10 12
table(enrichRCTM_top10$Cluster)
# C2 C5
#  5  5


# Ensure ONTOLOGY is factor with fixed order
enrichRCTM_top10$ONTOLOGY <- "Reactome"

# Plot
ggplot(enrichRCTM_top10, aes(x = FoldEnrichment,
                      y = reorder(Description, FoldEnrichment),
                      size = negLogFDR,
                      color = Cluster)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_size_continuous(name = expression(-log[10]~FDR)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  xlab("Fold enrichment") + ylab("Pathways")


##############################
## Wiki Pathways (enrichWP) ##
##############################

# Run enrichPathway for each cluster
enrich_results_WP <- lapply(unique(nom_annot_split$Cluster), function(cl){
  genes <- na.omit(nom_annot_split$EntrezGeneID[nom_annot_split$Cluster == cl])
  eg <- enrichWP(gene = genes,
                 organism = "Homo sapiens",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr"
                 )
  if(!is.null(eg)){
    as.data.frame(eg) %>% mutate(Cluster = cl)
  }
})
# Combine all clusters
enrich_WP <- bind_rows(enrich_results_WP)
table(enrich_WP$Cluster)
# C5
# 35

# Add enrichment score (fold enrichment) and significance
enrich_WP <- enrich_WP %>%
  mutate(negLogFDR = -log10(p.adjust))

enrichWP_top10 <- enrich_WP %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup()
dim(enrichWP_top10) # 10 12
table(enrichWP_top10$Cluster)
# C5
#  5


# Ensure ONTOLOGY is factor with fixed order
enrichWP_top10$ONTOLOGY <- "Wiki Pathways"

# Plot
ggplot(enrichWP_top10, aes(x = FoldEnrichment,
                      y = reorder(Description, FoldEnrichment),
                      size = negLogFDR,
                      color = Cluster)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_size_continuous(name = expression(-log[10]~FDR)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  xlab("Fold enrichment") + ylab("Pathways")


#################################
## Disease Ontology (enrichDO) ##
#################################

# Run enrichDO for each cluster
enrich_results_DO <- lapply(unique(nom_annot_split$Cluster), function(cl){
  genes <- na.omit(nom_annot_split$EntrezGeneID[nom_annot_split$Cluster == cl])
  eg <- enrichDO(gene = genes,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr",
                 readable = TRUE)
  if(!is.null(eg)){
    as.data.frame(eg) %>% mutate(Cluster = cl)
  }
})
# Combine all clusters
enrich_DO <- bind_rows(enrich_results_DO)
table(enrich_DO$Cluster)
# C2 C5
#  4 30

# Add enrichment score (fold enrichment) and significance
enrich_DO <- enrich_DO %>%
  mutate(negLogFDR = -log10(p.adjust))

enrichDO_top10 <- enrich_DO %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup()
dim(enrichDO_top10) # 9 14
table(enrichDO_top10$Cluster)
# C2 C5
#  4  5


# Ensure ONTOLOGY is factor with fixed order
enrichDO_top10$ONTOLOGY <- "Disease Ontology"

# Plot
ggplot(enrichDO_top10, aes(x = FoldEnrichment,
                      y = reorder(Description, FoldEnrichment),
                      size = negLogFDR,
                      color = Cluster)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_size_continuous(name = expression(-log[10]~FDR)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  xlab("Fold enrichment") + ylab("Pathways")
  

save(nom_annot_split, comp_annot_split, enrich_df, enrich_top10, enrich_KEGG, enrichKEGG_top10, 
  enrich_RCTM, enrichRCTM_top10, enrich_WP, enrichWP_top10, enrich_DO, enrichDO_top10, file="clusterProfiler_pathways.RData")

write.csv(enrich_df, file="enrich_GO.csv", row.names=F, quote=F)
write.csv(enrich_KEGG, file="enrich_KEGG.csv", row.names=F, quote=F)
write.csv(enrich_RCTM, file="enrich_Reactome.csv", row.names=F, quote=F)
write.csv(enrich_WP, file="enrich_WikiPathways.csv", row.names=F, quote=F)
write.csv(enrich_DO, file="enrich_DO.csv", row.names=F, quote=F)


#############################
## Top 5 Pathways Plotting ##
#############################

load("nom_comp_annot_split.RData")

extract_cols <- c("Description", "FoldEnrichment", "pvalue", "Cluster", "ONTOLOGY")

# GO
c1_GO <- enrichGO(gene = nom_annot_split[nom_annot_split$Cluster == "C1",]$EntrezGeneID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c1_GO@result$Cluster <- "C1"
c1_GO@result$ONTOLOGY <- paste0("GO:", c1_GO@result$ONTOLOGY)
dim(c1_GO) # 36 14
head(c1_GO@result[c1_GO@result$ONTOLOGY == "BP",], 5)
c1_GO_fig <- c1_GO@result[c1_GO@result$Description %in% c("response to fibroblast growth factor", "postsynapse organization", "cytokine activity"),][,extract_cols]

c2_GO <- enrichGO(gene = nom_annot_split[nom_annot_split$Cluster == "C2",]$EntrezGeneID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c2_GO@result$Cluster <- "C2"
c2_GO@result$ONTOLOGY <- paste0("GO:", c2_GO@result$ONTOLOGY)
dim(c2_GO) # 558  14
c2_GO_fig <- c2_GO@result[c2_GO@result$Description %in% c("leukocyte apoptotic process", "ubiquitin conjugating enzyme activity", "endopeptidase activity"),][,extract_cols]

c3_GO <- enrichGO(gene = nom_annot_split[nom_annot_split$Cluster == "C3",]$EntrezGeneID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c3_GO@result$Cluster <- "C3"
c3_GO@result$ONTOLOGY <- paste0("GO:", c3_GO@result$ONTOLOGY)
dim(c3_GO) # 86 14
c3_GO_fig <- c3_GO@result[c3_GO@result$Description %in% c("positive regulation of leukocyte activation", "positive regulation of T cell differentiation", "neuropeptide activity", "distal axon"),][,extract_cols]

c4_GO <- enrichGO(gene = nom_annot_split[nom_annot_split$Cluster == "C4",]$EntrezGeneID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c4_GO@result$Cluster <- "C4"
c4_GO@result$ONTOLOGY <- paste0("GO:", c4_GO@result$ONTOLOGY)
dim(c4_GO) # 143  14
c4_GO_fig <- c4_GO@result[c4_GO@result$Description %in% c("vesicle lumen", "plasma lipoprotein particle", "cytoplasmic stress granule", "cell chemotaxis", "sphingolipid biosynthetic process"),][,extract_cols]

c5_GO <- enrichGO(gene = nom_annot_split[nom_annot_split$Cluster == "C5",]$EntrezGeneID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c5_GO@result$Cluster <- "C5"
c5_GO@result$ONTOLOGY <- paste0("GO:", c5_GO@result$ONTOLOGY)
dim(c5_GO) # 813  14
c5_GO_fig <- c5_GO@result[c5_GO@result$Description %in% c("regulation of adaptive immune response", "peptidyl-tyrosine modification", "chemokine activity", "immune receptor activity", "neuronal cell body", "Golgi lumen", "integrin complex"),][,extract_cols]


# Reactome
c1_Reactome <- enrichPathway(gene = nom_annot_split[nom_annot_split$Cluster == "C1",]$EntrezGeneID, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c1_Reactome@result$Cluster <- "C1"
c1_Reactome@result$ONTOLOGY <- "Reactome"
c1_Reactome@result$FoldEnrichment <- parse_ratio(c1_Reactome@result$GeneRatio) / parse_ratio(c1_Reactome@result$BgRatio)
dim(c1_Reactome) # 8 12
c1_RCTM_fig <- c1_Reactome@result[c1_Reactome@result$Description %in% c("Signaling by Interleukins", "FGFR3 ligand binding and activation"),][,extract_cols]

c2_Reactome <- enrichPathway(gene = nom_annot_split[nom_annot_split$Cluster == "C2",]$EntrezGeneID, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c2_Reactome@result$Cluster <- "C2"
c2_Reactome@result$ONTOLOGY <- "Reactome"
c2_Reactome@result$FoldEnrichment <- parse_ratio(c2_Reactome@result$GeneRatio) / parse_ratio(c2_Reactome@result$BgRatio)
dim(c2_Reactome) # 89  12
c2_RCTM_fig <- c2_Reactome@result[c2_Reactome@result$Description %in% c("TCR signaling", "MTOR signalling"),][,extract_cols]

c3_Reactome <- enrichPathway(gene = nom_annot_split[nom_annot_split$Cluster == "C3",]$EntrezGeneID, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c3_Reactome@result$Cluster <- "C3"
c3_Reactome@result$ONTOLOGY <- "Reactome"
c3_Reactome@result$FoldEnrichment <- parse_ratio(c3_Reactome@result$GeneRatio) / parse_ratio(c3_Reactome@result$BgRatio)
dim(c3_Reactome) # 41  12
c3_RCTM_fig <- c3_Reactome@result[c3_Reactome@result$Description %in% c("Downregulation of TGF-beta receptor signaling", "p53-Dependent G1 DNA Damage Response"),][,extract_cols]

c4_Reactome <- enrichPathway(gene = nom_annot_split[nom_annot_split$Cluster == "C4",]$EntrezGeneID, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c4_Reactome@result$Cluster <- "C4"
c4_Reactome@result$ONTOLOGY <- "Reactome"
c4_Reactome@result$FoldEnrichment <- parse_ratio(c4_Reactome@result$GeneRatio) / parse_ratio(c4_Reactome@result$BgRatio)
dim(c4_Reactome) # 5 12
c4_RCTM_fig <- c4_Reactome@result[c4_Reactome@result$Description %in% c("Protein ubiquitination"),][,extract_cols]

c5_Reactome <- enrichPathway(gene = nom_annot_split[nom_annot_split$Cluster == "C5",]$EntrezGeneID, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c5_Reactome@result$Cluster <- "C5"
c5_Reactome@result$ONTOLOGY <- "Reactome"
c5_Reactome@result$FoldEnrichment <- parse_ratio(c5_Reactome@result$GeneRatio) / parse_ratio(c5_Reactome@result$BgRatio)
dim(c5_Reactome) # 184   12
c5_RCTM_fig <- c5_Reactome@result[c5_Reactome@result$Description %in% c("Signaling by TGFB family members", "Degradation of AXIN", "Autophagy"),][,extract_cols]


# WikiParhways
c1_WP <- enrichWP(gene = nom_annot_split[nom_annot_split$Cluster == "C1",]$EntrezGeneID, organism = "Homo sapiens", pvalueCutoff = 0.05, pAdjustMethod = "none")
c1_WP <- setReadable(c1_WP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
c1_WP@result$Cluster <- "C1"
c1_WP@result$ONTOLOGY <- "WikiPathways"
dim(c1_WP) # 17 14
c1_WP_fig <- c1_WP@result[c1_WP@result$Description %in% c("Cytokine cytokine receptor interaction", "Androgen biosynthesis"),][,extract_cols]

c2_WP <- enrichWP(gene = nom_annot_split[nom_annot_split$Cluster == "C2",]$EntrezGeneID, organism = "Homo sapiens", pvalueCutoff = 0.05, pAdjustMethod = "none")
c2_WP <- setReadable(c2_WP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
c2_WP@result$Cluster <- "C2"
c2_WP@result$ONTOLOGY <- "WikiPathways"
dim(c2_WP) # 20 14
c2_WP_fig <- c2_WP@result[c2_WP@result$Description %in% c("IL17 signaling", "Prostaglandin signaling"),][,extract_cols]

c5_WP <- enrichWP(gene = nom_annot_split[nom_annot_split$Cluster == "C5",]$EntrezGeneID, organism = "Homo sapiens", pvalueCutoff = 0.05, pAdjustMethod = "none")
c5_WP <- setReadable(c5_WP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
c5_WP@result$Cluster <- "C5"
c5_WP@result$ONTOLOGY <- "WikiPathways"
dim(c5_WP) # 100  14
c5_WP_fig <- c5_WP@result[c5_WP@result$Description %in% c("Cytokines and inflammatory response", "PI3K Akt signaling", "Autophagy"),][,extract_cols]


# KEGG
# Even without FDR adjustment, only C5 gives pathways
c5_KEGG <- enrichKEGG(gene = nom_annot_split[nom_annot_split$Cluster == "C5",]$EntrezGeneID, organism = "hsa", pvalueCutoff = 0.05, pAdjustMethod = "none", use_internal_data = FALSE)
c5_KEGG <- setReadable(c5_KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
c5_KEGG@result$Cluster <- "C5"
c5_KEGG@result$ONTOLOGY <- "KEGG"
dim(c5_KEGG) # 61 16
c5_KEG_fig <- c5_KEGG@result[c5_KEGG@result$Description %in% c("FoxO signaling pathway", "Longevity regulating pathway", "JAK-STAT signaling pathway", "Longevity regulating pathway - multiple species", "Cytokine-cytokine receptor interaction"),][,extract_cols]


# DO
# Only C2 and C5 give pathways for DO
c2_DO <- enrichDO(gene = nom_annot_split[nom_annot_split$Cluster == "C2",]$EntrezGeneID, pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c2_DO@result$Cluster <- "C2"
c2_DO@result$ONTOLOGY <- "DO"
dim(c2_DO) # 73 14
c2_DO_fig <- c2_DO@result[c2_DO@result$Description %in% c("Alzheimer's disease", "tauopathy", "neuropathy", "macular degeneration", "vascular dementia"),][,extract_cols]

c5_DO <- enrichDO(gene = nom_annot_split[nom_annot_split$Cluster == "C5",]$EntrezGeneID, pvalueCutoff = 0.05, pAdjustMethod = "none", readable = TRUE)
c5_DO@result$Cluster <- "C5"
c5_DO@result$ONTOLOGY <- "DO"
dim(c5_DO) # 154  14
c5_DO_fig <- c5_DO@result[c5_DO@result$Description %in% c("diabetic neuropathy", "blood protein disease", "hyperglycemia",  "central nervous system benign neoplasm", "cerebral infarction"),][,extract_cols]


ST_GO <- rbind(as.data.frame(c1_GO), as.data.frame(c2_GO), as.data.frame(c3_GO), as.data.frame(c4_GO), as.data.frame(c5_GO))
dim(ST_GO) # 1636   14

ST_Reactome <- rbind(as.data.frame(c1_Reactome), as.data.frame(c2_Reactome), as.data.frame(c3_Reactome), as.data.frame(c4_Reactome), as.data.frame(c5_Reactome))
dim(ST_Reactome) # 327  12

ST_WP <- rbind(as.data.frame(c1_WP), as.data.frame(c2_WP), as.data.frame(c5_WP))
dim(ST_WP) # 137  14

ST_KEGG <- as.data.frame(c5_KEGG)
dim(ST_KEGG) # 61 16

ST_DO <- rbind(as.data.frame(c2_DO), as.data.frame(c5_DO))
dim(ST_DO) # 227  14

ST_Pathways <- dplyr::bind_rows(ST_KEGG, ST_Reactome, ST_GO, ST_WP, ST_DO)
ST_Pathways$p.adjust <- NULL
dim(ST_Pathways) # 2388   15
ST_Pathways[is.na(ST_Pathways)] <- ""
write.table(ST_Pathways, file="SuppTable_Pathways_Nominal_10232025.txt", sep="\t", row.names=F, quote=F)

top_pathways <- rbind(c1_GO_fig, c2_GO_fig, c3_GO_fig, c4_GO_fig, c5_GO_fig,
  c1_RCTM_fig, c2_RCTM_fig, c3_RCTM_fig, c4_RCTM_fig, c5_RCTM_fig,
  c1_WP_fig, c2_WP_fig, c5_WP_fig, c5_KEG_fig, c2_DO_fig, c5_DO_fig)
dim(top_pathways) # 53  5
top_pathways$negLogP <-  -log10(top_pathways$pvalue)
top_pathways$ONTOLOGY <- ifelse(top_pathways$ONTOLOGY == "DO", "Disease Ontology", top_pathways$ONTOLOGY)
top_pathways$ONTOLOGY <- ifelse(top_pathways$ONTOLOGY == "WikiPathways", "Wiki Pathways", top_pathways$ONTOLOGY)
top_pathways$ONTOLOGY <- factor(top_pathways$ONTOLOGY, levels = c("GO:BP", "KEGG", "GO:CC", "Reactome", "GO:MF", "Wiki Pathways", "Disease Ontology"))

plot_top_pathways <- ggplot(top_pathways, aes(x = FoldEnrichment,
                      y = reorder(Description, FoldEnrichment),
                      size = negLogP,
                      color = Cluster)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 2) + # strip.position = "right" (if you want Ontology label on the right-side)
  scale_size_continuous(name = expression(-log[10]~P)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  xlab("Fold enrichment") + ylab("Pathways")

pdf("Top_5_pathways_dotPlot_v2.pdf", width = 12, height = 10)
print(plot_top_pathways)
dev.off()

png("Top_5_pathways_dotPlot_v2.png", units="mm", width=240, height=190, res=1000)
print(plot_top_pathways)
dev.off()

save(plot_top_pathways, top_pathways, file="Top_5_pathways_dotPlot_v2.RData")


#######################################################################################
## Pathway enrichment of 54 Concordant and Nominal proteins in K-ADRC and Bio-Hermes ##
#######################################################################################

# load pre-compile RData of 454 nominal and all 7K proteins sumstats that have been split on "|" to have only one protein per row
load("nom_comp_annot_split.RData")
dim(nom_annot_split) # 451   9
dim(comp_annot_split) # 6049    8

# Read 54 proteins nominal and concorant across K-ADRC and Bio-Hermes
nominal_54 <- read.csv("KADRC_BioHermes_Nominal_Concardant_ABpos_only_54proteins.csv")
dim(nominal_54)
head(nominal_54, 2)
table(nominal_54$Protein %in% nom_annot_split$Protein)
# TRUE
#   54
nominal_54 <- inner_join(nominal_54, nom_annot_split[,c("Protein", "UniProt", "EntrezGeneID", "EntrezGeneSymbol", "Cluster")], by="Protein")
dim(nominal_54) # 57 14
head(nominal_54, 2)
table(nominal_54$Cluster)
# C1 C2 C3 C4 C5
#  7  9  7  9 25

# Gene Ontology (EnrichGO)
enrichGO_GO <- enrichGO(gene = na.omit(nominal_54$EntrezGeneID),
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pAdjustMethod = "fdr",
                 readable = TRUE)
dim(enrichGO_GO) # 2 13
enrichGO_GO@result$ONTOLOGY <- "GO:CC"
GO_df <- enrichGO_GO@result[enrichGO_GO@result$p.adjust < 0.05,]
GO_df <- GO_df[,c("ID", "Description", "p.adjust", "geneID")]
GO_df$Database <- "GO"

enrich_KEGG <- enrichKEGG(gene = na.omit(nominal_54$UniProt),
                 organism = "hsa",
                 keyType = "uniprot",
                 pAdjustMethod = "fdr")
enrich_KEGG <- setReadable(enrich_KEGG, OrgDb = org.Hs.eg.db, keyType="UNIPROT")
dim(enrich_KEGG) # 19 14
enrich_KEGG@result$ONTOLOGY <- "KEGG"
head(enrich_KEGG)
KEGG_df <- enrich_KEGG@result[enrich_KEGG@result$p.adjust < 0.05,]
KEGG_df <- KEGG_df[,c("ID", "Description", "p.adjust", "geneID")]
KEGG_df$Database <- "KEGG"

enrich_REACTOME <- enrichPathway(gene = na.omit(nominal_54$EntrezGeneID),
                 organism = "human",
                 pAdjustMethod = "fdr")
enrich_REACTOME <- setReadable(enrich_REACTOME, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dim(enrich_REACTOME) # 28  9
enrich_REACTOME@result$ONTOLOGY <- "REACTOME"
head(enrich_REACTOME)
REACTOME_df <- enrich_REACTOME@result[enrich_REACTOME@result$p.adjust < 0.05,]
REACTOME_df <- REACTOME_df[,c("ID", "Description", "p.adjust", "geneID")]
REACTOME_df$Database <- "REACTOME"

enrich_WikiPathways <- enrichWP(gene = na.omit(nominal_54$EntrezGeneID),
                 organism = "Homo sapiens",
                 pAdjustMethod = "fdr")
enrich_WikiPathways <- setReadable(enrich_WikiPathways, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dim(enrich_WikiPathways) # 34 12
enrich_WikiPathways@result$ONTOLOGY <- "WikiPathways"
head(enrich_WikiPathways)
WikiPath_df <- enrich_WikiPathways@result[enrich_WikiPathways@result$p.adjust < 0.05,]
WikiPath_df <- WikiPath_df[,c("ID", "Description", "p.adjust", "geneID")]
WikiPath_df$Database <- "WikiPathways"

enrich_DO <- enrichDO(gene = na.omit(nominal_54$EntrezGeneID),
                 pAdjustMethod = "fdr",
                 readable = TRUE)
dim(enrich_DO) # 2 12
enrich_DO@result$ONTOLOGY <- "DO"
ST_KEGG <- as.data.frame(enrich_KEGG)
dim(ST_KEGG) # 19 15
ST_Reactome <- as.data.frame(enrich_REACTOME)
dim(ST_Reactome) # 28 10
ST_GO <- as.data.frame(enrichGO_GO)
dim(ST_GO) # 2 13
ST_WP <- as.data.frame(enrich_WikiPathways)
dim(ST_WP) # 34 13

ST_Pathways <- dplyr::bind_rows(ST_KEGG, ST_Reactome, ST_GO, ST_WP)
dim(ST_Pathways) # 83 15
ST_Pathways[is.na(ST_Pathways)] <- ""
write.table(ST_Pathways, file="SuppTable_Pathways_54proteins_03042026.txt", sep="\t", row.names=F, quote=F)
save(enrich_KEGG, enrich_REACTOME, enrichGO_GO, enrich_WikiPathways, enrich_DO, file="SuppTable_Pathways_54proteins_03042026.RData")

GO_df <- GO_df[GO_df$Description == "protein kinase complex",]
KEGG_df <- KEGG_df[KEGG_df$Description %in% c("FoxO signaling pathway",
  "Longevity regulating pathway", "Insulin signaling pathway",
  "Oxytocin signaling pathway", "Cell cycle"),]
REACTOME_df <- REACTOME_df[REACTOME_df$Description %in% c("Regulation of TP53 Activity",
  "Interferon Signaling", "Signaling by FGFR", "MTOR signalling", "Macroautophagy"),]
WikiPath_df <- WikiPath_df[WikiPath_df$Description %in% c("Age related macular degeneration",
  "Synaptic signaling associated with autism spectrum disorder", "Autophagy", "Lipid metabolism pathway",
  "PI3K Akt signaling"),]
pathway_df <- rbind(GO_df, KEGG_df, REACTOME_df, WikiPath_df)
dim(pathway_df) # 16  5

nominal_54_df <- nominal_54 %>%
  separate_rows(Symbol, sep = "\\|")
nominal_54_df <- unique(nominal_54_df[,c("B_KADRC", "Symbol")])
dim(nominal_54_df) # 57  2

## 1. Split geneID column to long format
path_long <- pathway_df %>%
  tidyr::separate_rows(geneID, sep = "/") %>%   # one gene per row
  dplyr::rename(Symbol = geneID)
dim(path_long) # 73  5
## 2. Keep only genes that have B_KADRC values, join effect size
dat_plot <- path_long %>%
  dplyr::left_join(
    nominal_54_df %>% dplyr::select(Symbol, B_KADRC),
    by = "Symbol"
  )
dim(dat_plot) # 73  6
## 3. Order factors (optional, for nicer axes)
dat_plot <- dat_plot %>%
  mutate(
    Description = factor(Description,
                         levels = rev(unique(Description))),  # y order
    Symbol = factor(Symbol,
                    levels = unique(Symbol))                  # x order
  )
dim(dat_plot) # 73  6
## 4. Heatmap with one facet per Database
pathway_fig <- ggplot(dat_plot, aes(x = Symbol, y = Description)) +
  # tiles: blue–white–red, color‑blind friendly-ish
  geom_tile(aes(fill = B_KADRC), color = "white") +
  # points: database, using Okabe–Ito palette (color‑blind safe)
  geom_point(aes(color = Database), size = 2) +
  scale_fill_gradient2(
    low  = "#0072B2",   # blue
    mid  = "white",
    high = "#D55E00",   # reddish‑orange
    midpoint = 0,
    name = "foldChange"
  ) +
  scale_color_manual(
    values = c(
      "GO"           = "#0072B2",
      "KEGG"         = "Black",
      "REACTOME"     = "#009E73",
      "WikiPathways" = "#CC79A7"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),                  # remove grey grid lines
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title  = element_blank(),
    legend.position = "right"
  )
# png("proteomic_signature_pathways_dotPlot.png", units="mm", width=240, height=190, res=1000)
# print(pathway_fig)
# dev.off()
ggsave("proteomic_signature_pathways_dotPlot.pdf", plot = pathway_fig, width = 12, height = 4, dpi = 600)
ggsave("proteomic_signature_pathways_dotPlot.png", pathway_fig, width = 12, height = 4, dpi = 600)


save(pathway_fig, dat_plot, file="proteomic_signature_pathways.RData")

geneList <- nominal_54$B_KADRC
names(geneList) <- nominal_54$EntrezGeneID
heatplot(enrichGO_GO, foldChange=geneList, showCategory=5)
heatplot(enrich_KEGG, foldChange=geneList, showCategory=5)
heatplot(enrich_REACTOME, foldChange=geneList, showCategory=5)

de <- names(geneList)[abs(geneList) > 2]
edo <- enrichDGN(de)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
plot_list(p1, p2, ncol=1, tag_levels = 'A')