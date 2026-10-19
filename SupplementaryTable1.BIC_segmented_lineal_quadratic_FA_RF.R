rm(list = ls())

library(openxlsx)
library(readxl)
library(dplyr)
library(readr)
library(rsample)
library(modelr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)

brain_size <- read.csv("D:/Ageing/Ampliación_paper_reviews/WM_MRI/brain_sizes.csv")

#FA
ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/FA/NoCleaned_FA_without_outliers/wm_data_FA_without_outliers.xlsx"
#FR
ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/FR/WM/wm_data_FR_without_outliers.xlsx"

gm_data <- read_excel(ruta_archivo)

gm_data$pnd <- as.numeric(gm_data$pnd)

# Los valores NA de pnd corresponden a 731 días
gm_data <- gm_data %>%
  mutate(
    pnd = if_else(is.na(pnd), 731, pnd)
  )

gm_data$region <- as.factor(gm_data$region)
gm_data$rat_code <- as.factor(gm_data$rat_code)
gm_data$sex <- as.factor(gm_data$sex)

# ============================================================
# 2. CENTRAR EDAD (Global para el modelo global)
# ============================================================
gm_data <- gm_data %>%
  mutate(
    pnd_c = pnd - mean(pnd, na.rm = TRUE)
  )
library(dplyr)

gm_data <- gm_data %>%
  left_join(
    brain_size,
    by = c("experiment" = "rat_code")
  )

library(dplyr)
library(segmented)

# ============================================================
# REGIONES
# ============================================================

gm_regiones <- unique(gm_data$region)

# ============================================================
# LISTA PARA GUARDAR RESULTADOS
# ============================================================

resultados_BIC <- list()

# ============================================================
# LOOP POR REGIÓN
# ============================================================

for (i in seq_along(gm_regiones)) {
  
  reg <- gm_regiones[i]
  
  cat("\nAnalizando:", reg, "\n")
  
  # ----------------------------------------------------------
  # Datos de la región
  # ----------------------------------------------------------
  
  datos_region <- gm_data %>%
    filter(region == reg) %>%
    filter(
      !is.na(FR),
      !is.na(pnd),
      !is.na(sex)
    )
  
  # ----------------------------------------------------------
  # MODELO CUADRÁTICO
  # ----------------------------------------------------------
  
  mod_quad <- lmer(
    FR ~ pnd_c * sex +
      I(pnd_c^2) + brain_size_voxels,
    data = datos_region
  )
  
  # ----------------------------------------------------------
  # MODELO LINEAL
  # ----------------------------------------------------------
  
  mod_lin <- lm(
    FR ~ pnd * sex + brain_size_voxels,
    data = datos_region
  )
  
  
  # ----------------------------------------------------------
  # MODELO SEGMENTADO
  # ----------------------------------------------------------
  
  mod_pw <- lm(
    FR ~ pnd * sex + brain_size_voxels,
    data = datos_region
  )
  
  seg_pw <- tryCatch(
    segmented(
      mod_pw,
      seg.Z = ~pnd,
      psi = 200
    ),
    error = function(e) NULL
  )
  
  # ----------------------------------------------------------
  # BIC
  # ----------------------------------------------------------
  
  bic_lin <- BIC(mod_lin)
  bic_quad <- BIC(mod_quad)
  
  if (!is.null(seg_pw)) {
    
    bic_seg <- BIC(seg_pw)
    
    breakpoint <- seg_pw$psi[, "Est."]
    
  } else {
    
    bic_seg <- NA_real_
    breakpoint <- NA_real_
  }
  
  # ----------------------------------------------------------
  # MODELO GANADOR
  # ----------------------------------------------------------
  
  bics <- c(
    Cuadratico = bic_quad,
    Lineal = bic_lin,
    Segmentado = bic_seg
  )
  
  modelo_ganador <- names(which.min(bics))
  
  # ----------------------------------------------------------
  # GUARDAR
  # ----------------------------------------------------------
  
  resultados_BIC[[i]] <- data.frame(
    region = reg,
    BIC_cuadratico = bic_quad,
    BIC_lineal = bic_lin,
    BIC_segmentado = bic_seg,
    breakpoint = breakpoint,
    modelo = modelo_ganador
  )
}

# ============================================================
# UNIR RESULTADOS
# ============================================================

resultados_BIC <- bind_rows(resultados_BIC)

resultados_BIC

write.csv(
  resultados_BIC,
  "LMM_results_WM_SEGMENTED_LINEAL_CUADRATICO_BIC_.csv",
  row.names = FALSE
)
