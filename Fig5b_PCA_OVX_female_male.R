# ============================================================
# PCA NULISA - PND161
# OVX vs Female vs Male
# Detectability filter includes OVX samples
# ============================================================

rm(list = ls())
# ============================================================
# LOAD ORIGINAL PROTEIN PANEL
# ============================================================

proteins_keep <- readRDS(
  "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/proteins_keep_Female_Male.rds"
)

# Verify
cat("Number of proteins in original panel:", length(proteins_keep), "\n")

# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(tidyverse)
library(readxl)
library(ggplot2)

# ============================================================
# 2. LOAD DATA
# ============================================================

file <- "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/Alamar_results/data/NULISAseq_NPQ_PP2025_0236_GLIA-MRI.xlsx"

df <- read_excel(
  file,
  sheet = "RQ NPQ Values"
)

metadata <- read_excel(
  "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/muestras_ensayo.xlsx",
  sheet = 2
)


# ============================================================
# 3. CLEANING AND DEFINE GROUPS
# ============================================================

df_final <- df %>%
  dplyr::select(
    SampleName,
    Target,
    NPQ,
    LOD
  ) %>%
  left_join(
    metadata,
    by = "SampleName"
  ) %>%
  # Exclude the four predefined samples
  filter(
    !SampleName %in% c(
      "B280222",
      "A230522",
      "B131222",
      "A300123"
    )
  ) %>%
  mutate(
    NPQ = as.numeric(NPQ),
    LOD = as.numeric(LOD),
    Sex = as.factor(Sex),
    Timepoint = as.numeric(Timepoint),
    Protein = Target,
    Group = case_when(
      startsWith(AnimalID, "O") ~ "OVX",
      as.character(Sex) == "female" ~ "Female",
      as.character(Sex) == "male" ~ "Male",
      TRUE ~ NA_character_
    )
  )

# Check group assignment
df_final %>%
  distinct(SampleName, AnimalID, Sex, Timepoint, Group) %>%
  count(Group)

# ============================================================
# 6. SELECT PND161
# ============================================================

pca_data <- df_final %>%
  filter(
    Timepoint == 161,
    Group %in% c("OVX", "Female", "Male"),
    Protein %in% proteins_keep
  )

# ============================================================
# 7. CHECK SAMPLE SIZE AT PND161
# ============================================================

sample_info <- pca_data %>%
  distinct(
    SampleName,
    AnimalID,
    Sex,
    Timepoint,
    Group
  ) %>%
  arrange(
    Group,
    AnimalID
  )

cat("\nSample counts at PND161:\n")

print(
  sample_info %>%
    count(Group)
)

cat("\nIndividual samples:\n")

print(sample_info)

# ============================================================
# 8. CREATE PROTEIN x SAMPLE MATRIX
# ============================================================

pca_wide <- pca_data %>%
  dplyr::select(
    SampleName,
    Protein,
    NPQ
  ) %>%
  pivot_wider(
    names_from = Protein,
    values_from = NPQ
  )

# Metadata for each sample
pca_metadata <- pca_data %>%
  dplyr::select(
    SampleName,
    AnimalID,
    Sex,
    Timepoint,
    Group
  ) %>%
  distinct() %>%
  filter(
    SampleName %in% pca_wide$SampleName
  ) %>%
  arrange(
    match(
      SampleName,
      pca_wide$SampleName
    )
  )

# Convert to numeric matrix
pca_matrix <- pca_wide %>%
  dplyr::select(-SampleName) %>%
  as.data.frame()

rownames(pca_matrix) <- pca_wide$SampleName

pca_matrix[] <- lapply(
  pca_matrix,
  as.numeric
)

# Check missing values
cat(
  "\nMissing values before removing proteins:",
  sum(is.na(pca_matrix)),
  "\n"
)

# Remove proteins with missing values
pca_matrix <- pca_matrix[
  ,
  colSums(is.na(pca_matrix)) == 0,
  drop = FALSE
]

# Identify proteins with zero or undefined variance
protein_sd <- sapply(
  pca_matrix,
  sd
)

zero_var <- !is.finite(protein_sd) | protein_sd == 0

cat(
  "\nProteins removed due to zero/undefined variance:",
  sum(zero_var),
  "\n"
)

print(
  names(pca_matrix)[zero_var]
)

# Remove constant or invalid proteins
pca_matrix <- pca_matrix[
  ,
  !zero_var,
  drop = FALSE
]

# Number of proteins entering PCA
cat(
  "\nNumber of proteins entering PCA:",
  ncol(pca_matrix),
  "\n"
)

if (ncol(pca_matrix) < 2) {
  stop("Fewer than two proteins remain for PCA.")
}

if (nrow(pca_matrix) < 2) {
  stop("Fewer than two samples remain for PCA.")
}

# ============================================================
# 9. PCA
# ============================================================

pca_result <- prcomp(
  pca_matrix,
  center = TRUE,
  scale. = TRUE
)

# Proteins excluded from PCA after detectability filtering
excluded_proteins <- setdiff(
  proteins_keep,
  colnames(pca_matrix)
)

cat(
  "\nProteins retained after detectability filtering:",
  length(proteins_keep),
  "\n"
)

cat(
  "Proteins excluded before PCA due to missing values,",
  "zero variance or undefined variance:",
  length(excluded_proteins),
  "\n"
)

print(excluded_proteins)

# Detailed information about excluded proteins
excluded_protein_data <- pca_data %>%
  filter(
    Protein %in% excluded_proteins
  ) %>%
  dplyr::select(
    SampleName,
    AnimalID,
    Group,
    Protein,
    NPQ,
    LOD
  ) %>%
  arrange(
    Protein,
    Group,
    AnimalID
  )

print(excluded_protein_data)

# ============================================================
# 10. VARIANCE EXPLAINED
# ============================================================

variance_explained <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

PC1_variance <- variance_explained[1]
PC2_variance <- variance_explained[2]

PC1_PC2_variance <- sum(
  variance_explained[1:2]
)

cat(
  "\nPC1 explains:",
  round(PC1_variance, 2),
  "% of total variance\n"
)

cat(
  "PC2 explains:",
  round(PC2_variance, 2),
  "% of total variance\n"
)

cat(
  "PC1 + PC2 explain:",
  round(PC1_PC2_variance, 2),
  "% of total variance\n"
)

# Full variance table
variance_table <- data.frame(
  PC = paste0(
    "PC",
    seq_along(variance_explained)
  ),
  Variance_explained = variance_explained
)

print(variance_table)

# ============================================================
# 11. PCA SCORES AND METADATA
# ============================================================

pca_scores <- as.data.frame(
  pca_result$x
)

pca_scores$SampleName <- rownames(pca_scores)

pca_scores <- pca_scores %>%
  left_join(
    pca_metadata,
    by = "SampleName"
  )

cat("\nPCA score preview:\n")

print(head(pca_scores))

cat("\nSamples represented in PCA:\n")

print(
  pca_scores %>%
    count(Group)
)

# ============================================================
# 12. PCA PLOT: PC1 vs PC2
# ============================================================

pca_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = Group
  )
) +
  geom_point(
    size = 3.5,
    alpha = 0.9
  ) +
  theme_classic(
    base_size = 15
  ) +
  labs(
    title = "NULISA PCA - PND161",
    x = paste0(
      "PC1 (",
      round(PC1_variance, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(PC2_variance, 1),
      "%)"
    ),
    color = "Group"
  )

print(pca_plot)

# ============================================================
# 13. PCA PLOT WITH 95% CONFIDENCE ELLIPSES
# ============================================================

pca_plot_ellipses <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    fill = Group
  )
) +
  stat_ellipse(
    aes(group = Group),
    type = "norm",
    level = 0.95,
    geom = "polygon",
    alpha = 0.10,
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.5,
    alpha = 0.9
  ) +
  theme_classic(
    base_size = 15
  ) +
  labs(
    title = "NULISA PCA - PND161",
    x = paste0(
      "PC1 (",
      round(PC1_variance, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(PC2_variance, 1),
      "%)"
    ),
    color = "Group",
    fill = "Group"
  )

print(pca_plot_ellipses)

