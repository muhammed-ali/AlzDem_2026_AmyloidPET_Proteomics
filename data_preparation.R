rm(list = ls())
options(stringsAsFactors = FALSE)

library(readr)
library(readxl)
library(data.table)
library(dplyr)
library(fuzzyjoin)
library(lubridate)
library(writexl)

#========================
# 1. Load SomaScan 7K plasma data
#========================

Zscore_Soma7k_Plasma <- read_csv("Zscore_Soma7k_Plasma_PIGEON_HAGERMAN_NIAGARA_COVID_Stanford_protein_matrix.csv")
MAP_Zscore <- Zscore_Soma7k_Plasma %>%
  dplyr::filter(grepl("^MAP", UniquePhenoID))

#========================
# 2. Load amyloid PET data
#========================

Amyloid <- read_excel("single_amyloid_for_each_visit_11_19_new 3.xlsx") %>%
  dplyr::mutate(
    UniquePhenoID = paste0("MAP_", ID),
    PET_Date = as.Date(PET_Date)
  )

#========================
# 3. Remove duplicated Soma samples by project priority
# Priority: Pigeon > Niagara > COVID
#========================

MAP_Zscore <- MAP_Zscore %>%
  dplyr::mutate(
    DrawDate = as.Date(DrawDate),
    Project_Priority = case_when(
      grepl("Pigeon", Project, ignore.case = TRUE) ~ 1,
      grepl("Niagara", Project, ignore.case = TRUE) ~ 2,
      grepl("COVID", Project, ignore.case = TRUE) ~ 3,
      TRUE ~ 4
    )
  ) %>%
  arrange(UniquePhenoID, DrawDate, Project_Priority) %>%
  group_by(UniquePhenoID, DrawDate) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(-Project_Priority)

#========================
# 4. Match Soma draw date with PET date within ±183 days
#========================

matched_all <- fuzzy_inner_join(
  MAP_Zscore %>% dplyr::select(UniquePhenoID, DrawDate),
  Amyloid %>% dplyr::select(UniquePhenoID, PET_Date),
  by = c(
    "UniquePhenoID" = "UniquePhenoID",
    "DrawDate" = "PET_Date"
  ),
  match_fun = list(
    `==`,
    function(x, y) abs(as.numeric(difftime(x, y, units = "days"))) <= 183
  )
) %>%
  dplyr::mutate(
    DayDiff = abs(as.numeric(difftime(DrawDate, PET_Date, units = "days")))
  )

# Keep the closest PET date for each Soma sample
matched_closest <- matched_all %>%
  group_by(UniquePhenoID.x, DrawDate) %>%
  slice_min(DayDiff, with_ties = FALSE) %>%
  ungroup()

#========================
# 5. Merge PET and Soma data
#========================

merged_part1 <- matched_closest %>%
  left_join(
    Amyloid,
    by = c("UniquePhenoID.y" = "UniquePhenoID", "PET_Date" = "PET_Date")
  )

merged_final <- merged_part1 %>%
  left_join(
    MAP_Zscore,
    by = c("UniquePhenoID.x" = "UniquePhenoID", "DrawDate" = "DrawDate")
  ) %>%
  dplyr::select(-UniquePhenoID.y) %>%
  dplyr::rename(UniquePhenoID = UniquePhenoID.x)

#========================
# 6. Create continuous amyloid z-score
# Priority: Centiloid_PiB > Centiloid_av45 > fsuvr_PiB > fsuvr_av45
#========================

continuous_pet_data <- merged_final %>%
  dplyr::mutate(
    combined_amyloid_zscore = coalesce(
      Centiloid_PiB_zscore,
      Centiloid_av45_zscore,
      fsuvr_PiB_zscore,
      fsuvr_av45_zscore
    ),
    tracer = case_when(
      !is.na(Centiloid_PiB_zscore) ~ "Centiloid_PiB",
      !is.na(Centiloid_av45_zscore) ~ "Centiloid_av45",
      !is.na(fsuvr_PiB_zscore) ~ "fsuvr_PiB",
      !is.na(fsuvr_av45_zscore) ~ "fsuvr_av45",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(combined_amyloid_zscore))

#========================
# 7. Select baseline sample per participant
#========================

continuous_pet_baseline <- continuous_pet_data %>%
  group_by(UniquePhenoID) %>%
  arrange(DrawDate) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(where(~ !all(is.na(.))))


#========================
# 8. Apply 85% call rate filter to proteins
#========================

protein_start <- which(grepl("^X", names(continuous_pet_baseline)))[1]

meta_cols <- continuous_pet_baseline[, 1:(protein_start - 1)]
protein_cols <- continuous_pet_baseline[, protein_start:ncol(continuous_pet_baseline)]

call_rate <- colMeans(!is.na(protein_cols))

failed_proteins <- names(call_rate[call_rate < 0.85])
passed_proteins <- protein_cols[, call_rate >= 0.85]

filtered_continuous_pet_baseline <- bind_cols(meta_cols, passed_proteins)

cat("Number of failed proteins:", length(failed_proteins), "\n")
cat("Final dimension after 85% call rate filter:\n")
print(dim(filtered_continuous_pet_baseline))

#========================
# 9. Calculate protein PCs
# Missing values are imputed by random sampling within each protein
#========================

prot <- filtered_continuous_pet_baseline %>%
  select(starts_with("X")) %>%
  as.data.frame()

row.names(prot) <- filtered_continuous_pet_baseline$ExtIdentifier

count_na <- apply(is.na(prot), 2, sum)

set.seed(2)

for (i in which(count_na != 0)) {
  index <- is.na(prot[, i])
  prot[index, i] <- sample(prot[!index, i], sum(index), replace = TRUE)
}

stopifnot(!any(is.na(prot)))

out_pc <- prcomp(prot, center = FALSE, scale = FALSE)

PC_scores <- as.data.frame(out_pc$x) %>%
  select(PC1, PC2, PC3) %>%
  mutate(ExtIdentifier = row.names(.))

#========================
# 10. Merge PCs and phenotype data
#========================

data_with_pcs <- merge(
  PC_scores,
  filtered_continuous_pet_baseline,
  by = "ExtIdentifier"
)

pheno <- read_excel("Phenotype_CSF_Plasma_09162025.xlsx") %>%
  dplyr::select(
    -Project, -ShortPhenoID, -DOB, -BirthYR,-PET_Date) %>%
  distinct(UniquePhenoID, DrawDate, .keep_all = TRUE)

final_data <- left_join(
  data_with_pcs,
  pheno,
  by = c("UniquePhenoID", "DrawDate")
)

# Move protein columns to the end
non_X_cols <- names(final_data)[!grepl("^X", names(final_data))]
X_cols <- names(final_data)[grepl("^X", names(final_data))]

final_data <- final_data[, c(non_X_cols, X_cols)]

#========================
# 11. Create final cleaned dataset for DAA
#========================

final_data <- final_data %>%
  mutate(
    Tracer = case_when(
      grepl("PiB", tracer, ignore.case = TRUE) ~ "PiB",
      grepl("av45", tracer, ignore.case = TRUE) ~ "av45",
      TRUE ~ NA_character_
    )
  )

subset_data <- final_data %>%
  select(
    UniquePhenoID, DrawDate, PET_Date, DayDiff,
    Tracer, tracer,
    combined_amyloid_zscore, final_decision,
    Age_at_draw, Sex, PC1, PC2,
    all_of(X_cols)
  )

write_xlsx(
  subset_data,
  "cleaned_data_558samples_6873proteins.xlsx"
)

cat("Final cleaned dataset dimension:\n")
print(dim(subset_data)) # 558 6885