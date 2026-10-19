###############################################################
# PROTEOMICS TRAJECTORIES - FINAL (MODEL CI VERSION)
###############################################################

rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(segmented)

# ============================================================
# DATA
# ============================================================

file <- "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/Alamar_results/data/NULISAseq_NPQ_PP2025_0236_GLIA-MRI.xlsx"

df <- read_excel(file, sheet = "RQ NPQ Values")
metadata <- read_excel("D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/muestras_ensayo.xlsx", sheet = 2)

df_final <- df %>%
  dplyr::select(SampleName, Target, NPQ, LOD) %>%
  left_join(metadata, by = "SampleName") %>%
  filter(
    !startsWith(AnimalID, "O"),
    !SampleName %in% c("B280222","A230522","B131222","A300123")
  ) %>%
  mutate(
    NPQ = as.numeric(NPQ),
    LOD = as.numeric(LOD),
    detected = NPQ > LOD,
    Sex = as.factor(Sex),
    Timepoint = as.numeric(Timepoint),
    AnimalID = as.factor(AnimalID),
    Protein = Target
  )

# ============================================================
# FILTER (LOD)
# ============================================================

min_group <- df_final %>%
  group_by(Sex) %>%
  summarise(n = n_distinct(SampleName), .groups="drop") %>%
  summarise(min(n)) %>%
  pull()

proteins_keep <- df_final %>%
  group_by(Protein) %>%
  summarise(det = sum(detected, na.rm=TRUE), .groups="drop") %>%
  filter(det >= min_group) %>%
  pull(Protein)

df_final <- df_final %>% filter(Protein %in% proteins_keep)

proteins <- unique(df_final$Protein)

# ============================================================
# OUTPUT
# ============================================================

results_final <- data.frame()

pdf("proteomics_trajectories_ALL_CI.pdf", width = 10, height = 6)

# ============================================================
# LOOP
# ============================================================

for (p in proteins) {
  
  cat("Processing:", p, "\n")
  
  data_p <- df_final %>%
    filter(Protein == p) %>%
    arrange(Timepoint)
  
  if(nrow(data_p) < 10 || length(unique(data_p$AnimalID)) < 2) next
  
  # ========================================================
  # MODELS
  # ========================================================
  
  # ------------------------
  # Linear model
  # ------------------------
  m_lin <- tryCatch(
    lmer(
      NPQ ~ Timepoint * Sex +
        (1 | AnimalID),
      data = data_p,
      REML = FALSE
    ),
    error = function(e) NULL
  )
  
  
  # ------------------------
  # Quadratic model
  # Timepoint * Sex
  # + Timepoint^2
  # ------------------------
  m_quad <- tryCatch(
    lmer(
      NPQ ~ Timepoint * Sex +
        I(Timepoint^2) +
        (1 | AnimalID),
      data = data_p,
      REML = FALSE
    ),
    error = function(e) NULL
  )
  
  # ========================================================
  # BIC MODEL SELECTION
  # ========================================================
  models <- list(
    Linear = m_lin,
    Quadratic = m_quad
  )
  
  models <- models[!vapply(models, is.null, logical(1))]
  
  if (length(models) == 0) next
  bic <- sapply(models, BIC)
  
  best <- names(which.min(bic))
  model <- models[[best]]
  
  
  # ========================================================
  # PVALUES
  # ========================================================
  
  a <- anova(model)
  
  age_p <- sex_p <- int_p <- NA
  
  
  # ------------------------
  # Linear
  # ------------------------
  
  if(best == "Linear"){
    
    age_p <- a["Timepoint", "Pr(>F)"]
    sex_p <- a["Sex", "Pr(>F)"]
    int_p <- a["Timepoint:Sex", "Pr(>F)"]
    
  }
  
  
  # ------------------------
  # Quadratic
  # ------------------------
  
  if(best == "Quadratic"){
    
    # Age = linear + quadratic effect
    age_p <- min(
      a[c("Timepoint", "I(Timepoint^2)"), "Pr(>F)"],
      na.rm = TRUE
    )
    
    sex_p <- a["Sex", "Pr(>F)"]
    
    # Only Timepoint × Sex interaction
    int_p <- a["Timepoint:Sex", "Pr(>F)"]
    
  }
  
  
  # ========================================================
  # EFFECT CLASSIFICATION
  # ========================================================
  
  effect <- case_when(
    
    age_p < 0.05 & sex_p < 0.05 & int_p < 0.05 ~
      "Age+Sex+Int",
    
    int_p < 0.05 ~
      "Interaction",
    
    age_p < 0.05 & sex_p < 0.05 ~
      "Age+Sex",
    
    age_p < 0.05 ~
      "Age",
    
    sex_p < 0.05 ~
      "Sex",
    
    TRUE ~
      "n.s."
  )
  
  
  # ========================================================
  # GRID + PREDICTION (MODEL CI)
  # ========================================================
  
  grid <- expand.grid(
    Timepoint = seq(
      min(data_p$Timepoint),
      max(data_p$Timepoint),
      length.out = 100
    ),
    Sex = levels(data_p$Sex)
  )
  
  
  pred <- predict(
    model,
    newdata = grid,
    se.fit = TRUE,
    re.form = NA
  )
  
  
  grid$pred <- pred$fit
  grid$upper <- pred$fit + 1.96 * pred$se.fit
  grid$lower <- pred$fit - 1.96 * pred$se.fit
  
  
  # ========================================================
  # SAVE TABLE
  # ========================================================
  
  results_final <- rbind(
    results_final,
    data.frame(
      Protein = p,
      Best_model = best,
      
      BIC_linear =
        ifelse(
          "Linear" %in% names(bic),
          bic["Linear"],
          NA
        ),
      
      BIC_quad =
        ifelse(
          "Quadratic" %in% names(bic),
          bic["Quadratic"],
          NA
        ),
      
      BIC_seg =
        ifelse(
          "Segmented" %in% names(bic),
          bic["Segmented"],
          NA
        ),
      
      Age_p = age_p,
      Sex_p = sex_p,
      Int_p = int_p,
      Effect = effect
    )
  )
  
  
  # ========================================================
  # PLOT (CI MODEL)
  # ========================================================
  
  print(
    ggplot(
      data_p,
      aes(Timepoint, NPQ, color = Sex)
    ) +
      
      geom_point(alpha = 0.5) +
      
      geom_ribbon(
        data = grid,
        aes(
          x = Timepoint,
          ymin = lower,
          ymax = upper,
          fill = Sex
        ),
        alpha = 0.15,
        inherit.aes = FALSE
      ) +
      
      geom_line(
        data = grid,
        aes(
          y = pred,
          group = Sex
        ),
        linewidth = 1.2
      ) +
      
      theme_bw() +
      
      ggtitle(
        paste0(
          p,
          "\nBest: ", best,
          " | Age p=", signif(age_p, 3),
          " Sex p=", signif(sex_p, 3),
          " Int p=", signif(int_p, 3)
        )
      )
  )
}

dev.off()

# ============================================================
# EXPORT
# ============================================================

write.csv(results_final,
          "Proteomics_model_results.csv",
          row.names = FALSE)

View(results_final)

# # =========================
# # SOLO NEFL
# # =========================
# ========================================================
# QUADRATIC MODEL — NEFL
# ========================================================
###AÑADIR EFECTO BATCH: 
library(dplyr)
library(lme4)
library(lmerTest)

data_nefl <- df_final %>%
  filter(
    Protein == "Nefl",
    !is.na(NPQ),
    !is.na(Timepoint),
    !is.na(Sex),
    !is.na(AnimalID)
  ) %>%
  mutate(
    Timepoint = as.numeric(as.character(Timepoint)),
    Batch = if_else(
      grepl("^[ES]", as.character(AnimalID)),
      "Batch2",
      "Batch1"
    ),
    Batch = factor(Batch)
  )


data_nefl$Timepoint <- as.numeric(as.character(data_nefl$Timepoint))


m_nefl_quad <- lmer(
  NPQ ~ Timepoint * Sex +
    I(Timepoint^2) + 
    (1 | AnimalID),
  data = data_nefl,
  REML = TRUE
)

summary(m_nefl_quad)
anova(m_nefl_quad)

library(emmeans)

emm_sex_time <- emmeans(
  m_nefl_quad,
  ~ Sex | Timepoint,
  at = list(Timepoint = c(161, 280, 364, 546, 609))
)

# Comparaciones entre sexos en cada edad
comp_sexo <- pairs(emm_sex_time, adjust = "BH")

# Estimaciones, IC 95 % y p ajustado
resultados <- as.data.frame(
  summary(comp_sexo, infer = c(TRUE, TRUE), level = 0.95)
)

resultados

library(ggplot2)

# 1. Generar predicciones continuas para dibujar curvas suaves
# Creamos una secuencia de tiempo desde el día mínimo (161) al máximo (609)
emm_curva <- emmeans(
  m_nefl_quad, 
  ~ Sex * Timepoint, 
  at = list(Timepoint = seq(161, 609, by = 10))
)

# 2. Convertir las predicciones del modelo a un Data Frame
df_pred <- as.data.frame(emm_curva)

# NOTA: emmeans genera automáticamente las columnas 'lower.CL' y 'upper.CL' (o 'asymp.LCL'/'asymp.UCL').
# Para evitar errores entre versiones de paquetes, recalculamos explícitamente el IC 95% aquí:
df_pred$upper_95 <- df_pred$emmean + 1.96 * df_pred$SE
df_pred$lower_95 <- df_pred$emmean - 1.96 * df_pred$SE

# 1. Generar predicciones continuas para dibujar curvas suaves
emm_curva <- emmeans(
  m_nefl_quad, 
  ~ Sex * Timepoint, 
  at = list(Timepoint = seq(161, 609, lenght.out = 100))
)

# 2. Convertir las predicciones del modelo a un Data Frame y calcular IC 95%
df_pred <- as.data.frame(emm_curva)
df_pred$upper_95 <- df_pred$emmean + 1.96 * df_pred$SE
df_pred$lower_95 <- df_pred$emmean - 1.96 * df_pred$SE

# 3. Graficar con ggplot2 ampliando los límites del eje X
p <- ggplot(df_pred, aes(x = Timepoint, y = emmean, color = Sex, fill = Sex)) +
  # Puntos reales alineados exactamente en su Timepoint
  geom_point(
    data = data_nefl, 
    aes(y = NPQ), 
    alpha = 0.4, 
    size = 2
  ) +
  # Sombreado con el Intervalo de Confianza al 95%
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.15, color = NA) +
  # Líneas curvas ajustadas por el modelo cuadrático
  geom_line(linewidth = 1.2) +
  # Paleta de colores personalizada con nombres en minúsculas
  scale_color_manual(values = c("male" = "#77DD77", "female" = "#9370DB"), labels = c("male" = "Male", "female" = "Female")) +
  scale_fill_manual(values = c("male" = "#77DD77", "female" = "#9370DB"), labels = c("male" = "Male", "female" = "Female")) +
  
  # MODIFICACIÓN CLAVE: Forzamos los límites y marcas exactas del eje X
  scale_x_continuous(
    breaks = c(160, 200, 300, 400, 500, 600), # Muestra solo tus días reales de muestreo
    limits = c(140, 630)                 # Da margen a la izquierda y derecha para que no se corte
  ) +
  
  # Etiquetas
  labs(
    title = "Evolución de Nefl (NPQ) a lo largo del tiempo",
    subtitle = "Modelo mixto cuadrático con datos individuales",
    x = "Timepoint (Días)",
    y = "NPQ",
    color = "Sexo",
    fill = "Sexo"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# Forzar el reinicio del dispositivo gráfico
dev.off()

# Mostrar el gráfico por pantalla
print(p)

