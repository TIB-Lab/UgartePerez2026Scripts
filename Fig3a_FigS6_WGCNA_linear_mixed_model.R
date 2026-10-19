
# ============================================================
# 1. INITIALIZATION
# ============================================================

rm(list = ls())

library(WGCNA)
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(pheatmap)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

dir.create("../Results/WGCNA", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. LOAD DATA
# ============================================================

file <- paste0(
  "D:/Ageing/BLOOD BIOMARKERS/NULISA/",
  "Nulisa_analysis/Alamar_results/data/",
  "NULISAseq_NPQ_PP2025_0236_GLIA-MRI.xlsx"
)

df <- read_excel(file, sheet = "RQ NPQ Values")

metadata <- read_excel(
  "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/muestras_ensayo.xlsx",
  sheet = 2
)
# ============================================================
# 3. CLEANING
# ============================================================

df_final <- df %>%
  dplyr::select(SampleName, Target, NPQ, LOD) %>%
  left_join(metadata, by = "SampleName") %>%
  filter(
    !startsWith(AnimalID, "O"),
    !SampleName %in% c(
      "B280222", "A230522", "B131222", "A300123"
    )
  ) %>%
  mutate(
    NPQ = as.numeric(NPQ),
    LOD = as.numeric(LOD),
    Sex = factor(Sex),
    Timepoint = as.numeric(Timepoint),
    Protein = Target
  )

# Comprobar que cada muestra tiene un solo registro de metadata
metadata_check <- df_final %>%
  distinct(SampleName, AnimalID, Sex, Timepoint) %>%
  count(SampleName) %>%
  filter(n > 1)

if (nrow(metadata_check) > 0) {
  stop("Hay muestras con metadata inconsistente. Revisa SampleName.")
}

# ============================================================
# 3.1. FILTER PROTEINS BY DETECTABILITY
# ============================================================

min_group_size <- df_final %>%
  distinct(SampleName, Sex) %>%
  count(Sex) %>%
  pull(n) %>%
  min()

df_final <- df_final %>%
  mutate(detected = NPQ > LOD)

proteins_keep <- df_final %>%
  group_by(Protein) %>%
  summarise(
    detected_samples = sum(detected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(detected_samples >= min_group_size) %>%
  pull(Protein)

df_final <- df_final %>%
  filter(Protein %in% proteins_keep)

# ============================================================
# 3.2. SAMPLE SIZE
# ============================================================

sample_info <- df_final %>%
  distinct(SampleName, AnimalID, Sex, Timepoint)

print(
  sample_info %>%
    summarise(
      n_samples = n_distinct(SampleName),
      n_animals = n_distinct(AnimalID)
    )
)

print(
  sample_info %>%
    group_by(Timepoint, Sex) %>%
    summarise(
      n_samples = n_distinct(SampleName),
      n_animals = n_distinct(AnimalID),
      .groups = "drop"
    ) %>%
    arrange(Timepoint, Sex)
)

# ============================================================
# 4. WIDE FORMAT
# IMPORTANT: retain AnimalID for the mixed models
# ============================================================

df_wide <- df_final %>%
  dplyr::select(
    SampleName, AnimalID, Protein, NPQ, Sex, Timepoint
  ) %>%
  pivot_wider(
    names_from = Protein,
    values_from = NPQ
  )

# ============================================================
# 5. EXPRESSION MATRIX
# ============================================================

datExpr <- df_wide %>%
  dplyr::select(
    -SampleName, -AnimalID, -Sex, -Timepoint
  ) %>%
  as.data.frame()

rownames(datExpr) <- df_wide$SampleName
datExpr[] <- lapply(datExpr, as.numeric)

# ============================================================
# 6. QUALITY CONTROL
# ============================================================

gsg <- goodSamplesGenes(datExpr, verbose = 3)
print(gsg$allOK)

if (!gsg$allOK) {
  datExpr <- datExpr[
    gsg$goodSamples,
    gsg$goodGenes,
    drop = FALSE
  ]
}

# ============================================================
# 7. SAMPLE CLUSTERING
# ============================================================

sampleTree <- hclust(dist(datExpr), method = "average")

plot(
  sampleTree,
  main = "Sample clustering",
  sub = "",
  xlab = ""
)

# ============================================================
# 8. SOFT THRESHOLD
# ============================================================

powers <- 1:20

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

# ============================================================
# 9. SCALE-FREE TOPOLOGY
# ============================================================

pdf("../Results/WGCNA/WGCNA_scale_free_topology.pdf",
    width = 6, height = 5)

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  type = "b",
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit",
  main = "Scale independence"
)

abline(h = 0.8, col = "red")
dev.off()

# ============================================================
# 10. MEAN CONNECTIVITY
# ============================================================

pdf("../Results/WGCNA/WGCNA_mean_connectivity.pdf",
    width = 6, height = 5)

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  type = "b",
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  main = "Mean connectivity"
)

dev.off()

# ============================================================
# 11. SELECT SOFT POWER
# ============================================================

softPower <- 9

# ============================================================
# 12. ADJACENCY MATRIX
# ============================================================

adjacency <- WGCNA::adjacency(
  datExpr,
  power = softPower,
  type = "signed"
)

k <- softConnectivity(datExpr, power = 6, type = "signed")
print(mean(k))

# ============================================================
# 13. TOM MATRIX
# ============================================================

TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# ============================================================
# 14. PROTEIN CLUSTERING
# ============================================================

geneTree <- hclust(as.dist(dissTOM), method = "average")

plot(
  geneTree,
  xlab = "",
  sub = "",
  main = "Protein clustering on TOM-based dissimilarity",
  labels = FALSE,
  hang = 0.04
)

# ============================================================
# 15. MODULE DETECTION
# ============================================================

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 8
)

print(table(dynamicMods))

# ============================================================
# 16. MODULE COLORS
# ============================================================

moduleColors <- labels2colors(dynamicMods)
print(table(moduleColors))

# ============================================================
# 17. DENDROGRAM + MODULES
# ============================================================

pdf("../Results/WGCNA/WGCNA_modules.pdf",
    width = 12, height = 6)

plotDendroAndColors(
  geneTree,
  moduleColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

dev.off()

# ============================================================
# 18. MODULE EIGENPROTEINS
# ============================================================

MEList <- moduleEigengenes(
  datExpr,
  colors = moduleColors
)

MEs <- MEList$eigengenes

# ============================================================
# 18.1. MERGE CLOSE MODULES
# ============================================================

MEDiss <- 1 - cor(MEs, use = "pairwise.complete.obs")
METree <- hclust(as.dist(MEDiss), method = "average")

plot(
  METree,
  main = "Clustering of module eigengenes",
  xlab = "",
  sub = ""
)

abline(h = 0.25, col = "red")

merge <- mergeCloseModules(
  datExpr,
  moduleColors,
  cutHeight = 0.25,
  verbose = 3
)

moduleColors <- merge$colors
MEs <- merge$newMEs

print(table(moduleColors))

pdf("../Results/WGCNA/WGCNA_merged_modules.pdf",
    width = 12, height = 6)

plotDendroAndColors(
  geneTree,
  moduleColors,
  "Merged module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

dev.off()

# ============================================================
# 19. TRAIT DATA FOR MIXED MODELS
# ============================================================

traitData <- df_wide %>%
  dplyr::select(SampleName, AnimalID, Sex, Timepoint) %>%
  distinct() %>%
  filter(SampleName %in% rownames(datExpr)) %>%
  arrange(match(SampleName, rownames(datExpr)))

stopifnot(all(traitData$SampleName == rownames(datExpr)))

traitData <- traitData %>%
  mutate(
    AnimalID = factor(AnimalID),
    Sex = droplevels(factor(Sex)),
    Timepoint = as.numeric(Timepoint)
  )

# Define Male as the reference level when the metadata uses
# either "male" or "Male".
sex_levels <- levels(traitData$Sex)
male_level <- sex_levels[tolower(sex_levels) == "male"]

if (length(male_level) == 1) {
  traitData$Sex <- relevel(traitData$Sex, ref = male_level)
}

if (nlevels(traitData$Sex) < 2) {
  stop("Sex tiene menos de dos niveles tras filtrar las muestras.")
}

if (sd(traitData$Timepoint, na.rm = TRUE) == 0) {
  stop("Timepoint no varía; no se puede estimar su efecto.")
}

# Standardize time to make its coefficient represent a change
# per one SD of time.
traitData$Timepoint_z <- as.numeric(scale(traitData$Timepoint))

# Align MEs exactly with the samples used in traitData
MEs <- MEs[traitData$SampleName, , drop = FALSE]

stopifnot(identical(rownames(MEs), traitData$SampleName))

# ============================================================
# 20. MODULE-TRAIT MIXED MODELS
# Outcome: standardized module eigengene
# Fixed effects: Sex + Timepoint
# Random intercept: AnimalID
# ============================================================

library(lme4)
library(lmerTest)
library(dplyr)
library(purrr)
library(tibble)

results_list <- list()

for (mod in colnames(MEs)) {
  
  # Incluye MEgrey, si está presente
  dat <- traitData %>%
    mutate(ME = as.numeric(MEs[, mod])) %>%
    select(AnimalID, Sex, Timepoint_z, ME) %>%
    filter(
      complete.cases(AnimalID, Sex, Timepoint_z, ME),
      is.finite(ME),
      is.finite(Timepoint_z)
    )
  
  # Comprobaciones mínimas
  if (nrow(dat) < 5 ||
      dplyr::n_distinct(dat$AnimalID) < 2 ||
      sd(dat$ME) == 0) {
    next
  }
  
  # Estandarizar el eigengene
  dat$ME_z <- as.numeric(scale(dat$ME))
  
  # Modelo longitudinal: efectos fijos de sexo y tiempo,
  # intercepto aleatorio por animal
  fit <- tryCatch(
    lmerTest::lmer(
      ME_z ~ Sex + Timepoint_z + (1 | AnimalID),
      data = dat,
      REML = TRUE
    ),
    error = function(e) {
      message("Model failed for ", mod, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(fit)) next
  
  coefs <- coef(summary(fit))
  terms_to_extract <- intersect(
    c("Sexfemale", "Timepoint_z"),
    rownames(coefs)
  )
  
  if (length(terms_to_extract) == 0) next
  
  # Intervalos de confianza del 95 % para los efectos fijos
  ci <- tryCatch(
    confint(
      fit,
      parm = terms_to_extract,
      method = "Wald"
    ),
    error = function(e) NULL
  )
  
  for (term in terms_to_extract) {
    
    estimate <- coefs[term, "Estimate"]
    se       <- coefs[term, "Std. Error"]
    df       <- coefs[term, "df"]
    t_value  <- coefs[term, "t value"]
    p_value  <- coefs[term, "Pr(>|t|)"]
    
    # IC del 95 %; NA si no se pudo calcular
    ci_low <- NA_real_
    ci_high <- NA_real_
    
    if (!is.null(ci) && term %in% rownames(ci)) {
      ci_low  <- ci[term, 1]
      ci_high <- ci[term, 2]
    }
    
    results_list[[length(results_list) + 1]] <- tibble(
      Module = mod,
      Term = term,
      Effect = estimate,
      SE = se,
      df = df,
      t_value = t_value,
      CI95_low = ci_low,
      CI95_high = ci_high,
      p_value = p_value
    )
  }
}

# Combinar resultados y corregir por comparaciones múltiples
effect_results <- bind_rows(results_list) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    Trait = case_when(
      Term == "Sexfemale" ~ "Sex (female vs male)",
      Term == "Timepoint_z" ~ "Timepoint (per SD)",
      TRUE ~ Term
    ),
    FDR_label = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(FDR)

print(effect_results)

# Exportar tabla completa para el material suplementario
write.csv(
  effect_results,
  "WGCNA_mixed_model_results.csv",
  row.names = FALSE
)
# ============================================================
# 21. HEATMAP: EFFECT SIZE + FDR
# ============================================================

# Keep the same symmetric color scale around zero
max_effect <- max(abs(effect_results$Effect), na.rm = TRUE)

if (!is.finite(max_effect) || max_effect == 0) {
  max_effect <- 1
}

effect_results$Module <- factor(
  effect_results$Module,
  levels = rev(unique(effect_results$Module))
)

effect_results$Trait <- factor(
  effect_results$Trait,
  levels = unique(effect_results$Trait)
)


library(dplyr)
library(ggplot2)

# Preparar resultados de los modelos mixtos
plot_data <- effect_results %>%
  mutate(
    Trait = recode(
      Trait,
      "SexFemale" = "Sex",
      "Sex" = "Sex",
      "Timepoint_z" = "Time / age"
    ),
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ "ns"
    ),
    Label = paste0(
      sprintf("%.2f", Effect),
      "\nFDR=", format.pval(FDR, digits = 2, eps = 0.001),
      "\n", Significance
    )
  )

# Heatmap de effect sizes
p_effect <- ggplot(
  plot_data,
  aes(x = Trait, y = Module, fill = Effect)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = 3) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Effect size"
  ) +
  labs(
    title = "Module-trait associations: mixed models",
    subtitle = "Efectos ajustados; asteriscos basados en FDR",
    x = NULL,
    y = "WGCNA module"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold")
  )

print(p_effect)

ggsave(
  "../Results/WGCNA/mixed_model_effect_sizes.png",
  p_effect,
  width = 7,
  height = max(5, 0.35 * n_distinct(plot_data$Module)),
  dpi = 300
)


ggsave(
  "../Results/WGCNA/WGCNA_module_trait_mixed_model_heatmap.pdf",
  plot = p_heatmap,
  width = 8,
  height = 6
)

ggsave(
  "../Results/WGCNA/WGCNA_module_trait_mixed_model_heatmap.png",
  plot = p_heatmap,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 22. IDENTIFY MODULES CONTAINING NPTX2 AND NEFL
# ============================================================

target_proteins <- c("Nptx2", "Nefl")

protein_modules_of_interest <- data.frame(
  Protein = target_proteins,
  Module = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(target_proteins)) {
  
  protein <- target_proteins[i]
  
  if (protein %in% colnames(datExpr)) {
    
    protein_module <- moduleColors[
      match(protein, colnames(datExpr))
    ]
    
    protein_modules_of_interest$Module[i] <- protein_module
    
  } else {
    message(protein, " is not present in datExpr.")
  }
}

print(protein_modules_of_interest)


# ============================================================
# 23. MODULE ASSIGNMENT FOR ALL PROTEINS
# ============================================================

module_assignment <- data.frame(
  Protein = colnames(datExpr),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

print(module_assignment)


# ============================================================
# 24. PLOT MODULE EIGENGENES FOR NPTX2 AND NEFL ONLY
# ============================================================

library(dplyr)
library(ggplot2)

proteins_to_plot <- c("Nptx2", "Nefl")

sample_metadata <- df_wide %>%
  dplyr::select(SampleName, Sex, Timepoint) %>%
  distinct()

for (protein in proteins_to_plot) {
  
  # Check whether the protein is present
  if (!protein %in% colnames(datExpr)) {
    message(protein, " is not present in datExpr.")
    next
  }
  
  # Identify the protein's module
  module <- moduleColors[
    match(protein, colnames(datExpr))
  ]
  
  ME_name <- paste0("ME", module)
  
  if (!ME_name %in% colnames(MEs)) {
    message("Module eigengene ", ME_name, " not found.")
    next
  }
  
  # Prepare plotting data
  ME_df <- data.frame(
    SampleName = rownames(MEs),
    ModuleEigengene = MEs[, ME_name]
  ) %>%
    left_join(sample_metadata, by = "SampleName") %>%
    filter(!is.na(Sex))
  
  # Plot
  p <- ggplot(
    ME_df,
    aes(x = Sex, y = ModuleEigengene, fill = Sex)
  ) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      alpha = 0.85,
      linewidth = 0.8
    ) +
    geom_jitter(
      width = 0.08,
      size = 1.8,
      color = "grey20",
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "female" = "#E9A3A3",
        "male" = "#8FC9B5"
      )
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank()
    ) +
    labs(
      title = paste0(protein, " module (", module, ")"),
      x = NULL,
      y = "Module eigengene"
    )
  
  print(p)
}


