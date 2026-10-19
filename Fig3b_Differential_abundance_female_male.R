###########ANALISIS DIFERENCIAL 
rm(list=ls())

library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)
library(limma)
library(ggplot2)

# -----------------------------
# 1. CARGA DATOS
# -----------------------------
file <- "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/Alamar_results/data/NULISAseq_NPQ_PP2025_0236_GLIA-MRI.xlsx"

df <- read_excel(
  file,
  sheet = "RQ NPQ Values"
)

metadata <- read_excel(
  "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/muestras_ensayo.xlsx",
  sheet = 2
)

# -----------------------------
# 2. LIMPIEZA Y METADATA
# -----------------------------
df_clean <- df %>%
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
  filter(
    !startsWith(AnimalID, "O"),
    !SampleName %in% c(
      "B280222",
      "A230522",
      "B131222",
      "A300123"
    ),
    Timepoint %in% c(161, 280, 364, 546, 609),
    Sex %in% c("male", "female")
  ) %>%
  mutate(
    NPQ = as.numeric(NPQ),
    LOD = as.numeric(LOD),
    Sex = factor(Sex),
    Timepoint = factor(Timepoint),
    Protein = Target
  )

# -----------------------------
# 3. FILTRO DE DETECTABILIDAD
# -----------------------------

min_group_size <- df_clean %>%
  distinct(SampleName, Sex) %>%
  count(Sex) %>%
  pull(n) %>%
  min()

df_clean <- df_clean %>%
  mutate(
    detected = NPQ > LOD
  )

proteins_keep <- df_clean %>%
  group_by(Protein) %>%
  summarise(
    detected_samples = sum(detected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    detected_samples >= min_group_size
  ) %>%
  pull(Protein)

df_clean <- df_clean %>%
  filter(
    Protein %in% proteins_keep
  )


# ============================================================
# SAVE ORIGINAL PROTEIN PANEL
# Female vs Male detectability filter
# ============================================================

saveRDS(
  proteins_keep,
  file = "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/proteins_keep_Female_Male.rds"
)

# Verify
length(proteins_keep)

# -----------------------------
# 3. DEFINIR OVX
# -----------------------------

df_filtrado <- df_clean
df_filtrado$Sex <- factor(df_filtrado$Sex)

# =========================================================
# 4. ANÁLISIS DIFERENCIAL POR EDAD
# =========================================================

df_filtrado <- df_clean %>%
  filter(Timepoint == 546) %>%   # cambiar a 280, 364 o 546
  mutate(
    Sex = factor(Sex)
  )

# -----------------------------
# 4. LONG FORMAT
# -----------------------------
# df_long <- df_filtrado %>%
#   pivot_longer(
#     cols = -c(SampleName, AnimalID, Timepoint, Sex, Observations),
#     names_to = "Protein",
#     values_to = "NPQ"
#   )

# -----------------------------
# 5. WIDE FORMAT (proteínas x muestras)
# -----------------------------
heatmap_data <- df_filtrado %>%
  dplyr::select(SampleName, Protein, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) 

mat <- as.data.frame(heatmap_data)
rownames(mat) <- mat$Protein
mat$Protein <- NULL

mat <- as.matrix(mat)
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- heatmap_data$Protein

# -----------------------------
# 6. LIMPIEZA 
# -----------------------------
mat <- mat[complete.cases(mat), ]
mat <- mat[apply(mat, 1, sd, na.rm = TRUE) != 0, ]

# =========================================================
# 📊 8. LIMMA 
# =========================================================
expr <- mat  

# Una fila de metadata por muestra y en el mismo orden que expr
sample_metadata <- df_filtrado %>%
  distinct(SampleName, AnimalID, Timepoint, Sex) %>%
  slice(match(colnames(expr), SampleName))

group <- factor(sample_metadata$Sex)

# Comprobaciones
ncol(expr)
length(group)
identical(colnames(expr), sample_metadata$SampleName)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts(
  Male_vs_Female = male - female,
  #OVX_vs_Female = OVX - female,
  #OVX_vs_Male = OVX - male,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res_mf <- topTable(fit2, coef="Male_vs_Female", number=Inf)
#res_ovx_f <- topTable(fit2, coef="OVX_vs_Female", number=Inf)
#res_ovx_m <- topTable(fit2, coef="OVX_vs_Male", number=Inf)

# -----------------------------
# 9. VOLCANO PLOT
# -----------------------------
# res_mf <- res_mf %>%
#   mutate(
#     Significant = adj.P.Val < 0.05,
#     neg_log10_FDR = -log10(adj.P.Val)
#   )
# 
# p <- ggplot(res_mf, aes(x = logFC, y = neg_log10_FDR)) +
#   geom_point(aes(color = Significant), size = 2) +
#   scale_color_manual(values = c("grey70", "red")) +
#   theme_classic(base_size = 14) +
#   xlab("log2 fold change (Male vs Female)") +
#   ylab("-log10(FDR)") +
#   ggtitle("Male vs Female")
# 
# ggsave("volcano_Male_vs_Female.png", p, width = 6, height = 5, dpi = 300)

res_mf$Protein <- rownames(res_mf)

res_mf <- res_mf %>%
  mutate(
    regulation = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Up in Male",
      adj.P.Val < 0.05 & logFC < 0 ~ "Up in Female",
      TRUE ~ "n.s."
    ),
    neg_log10_FDR = -log10(adj.P.Val)
  )

library(ggrepel)
ggplot(res_mf, aes(x = logFC, y = neg_log10_FDR)) +
  geom_point(
    aes(color = regulation),
    size = 2.5,
    alpha = 0.85
  ) +
  geom_text_repel(
    data = res_mf %>% filter(regulation != "n.s."),
    aes(label = Protein),
    size = 4,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Up in Male" = "firebrick",
      "Up in Female" = "royalblue",
      "n.s." = "grey75"
    ),
    name = "Differential abundance"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_classic(base_size = 14) +
  labs(
    title = "Male vs Female",
    x = "log2 fold change",
    y = "-log10(FDR)"
  )

ggsave(
  filename = "../Results/volcano_Male_vs_Female.pdf",
  width = 7,
  height = 6
)
write.csv(
  res_mf,
  "../Results/limma_Male_vs_Female_546PND.csv",
  row.names = FALSE
)



