rm(list=ls())
library(lme4)
library(lmerTest)
library(emmeans)

# Especificar la ruta del archivo xlsx
#MD
ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/MD/CleanedMD_without_Outliers/gm_data_MD_without_outliers.xlsx"

# Cargar el archivo xlsx en un data frame
wm_data <- read_excel(ruta_archivo)
wm_data$pnd <- as.numeric(wm_data$pnd)

#LOS VALORES DE 731 DÍAS PONE NA, SOLUCION:#
wm_data <- wm_data %>%
  mutate(pnd = if_else(is.na(pnd), 731, pnd))

wm_data$region <- as.factor(wm_data$region) #primero debemos factorizar la columna 'region'
wm_regiones <- levels(wm_data$region)

wm_data$sex_num <- as.numeric(wm_data$sex == 'Male')

str(wm_data)

table(wm_data$sex)
table(wm_data$region)
length(unique(wm_data$rat_code))

wm_data$rat_code <- as.factor(wm_data$rat_code)
wm_data$sex <- as.factor(wm_data$sex)

brain_size <- read.csv(
  "D:/Ageing/Ampliación_paper_reviews/WM_MRI/brain_sizes.csv"
)

# Unir brain size a wm_data
wm_data <- wm_data %>%
  left_join(brain_size, by = c("experiment" = "rat_code"))

resultados_BIC <- data.frame()

modelos_lineales <- list()
modelos_cuadraticos <- list()

for (reg in wm_regiones) {
  
  datos_region <- wm_data %>%
    filter(region == reg) %>%
    mutate(
      pnd_reg_c = pnd - mean(pnd, na.rm = TRUE)
    )
  
  # ==========================
  # Modelo lineal
  # ==========================
  
  mod_linear <- tryCatch(
    lmer(
      MD ~ pnd_reg_c * sex + brain_size_voxels +
        (1 | rat_code),
      data = datos_region,
      REML = FALSE
    ),
    error = function(e) NULL
  )
  
  
  # ==========================
  # Modelo cuadrático
  # ==========================
  
  mod_quad <- tryCatch(
    lmer(
      MD ~ pnd_reg_c * sex +
        I(pnd_reg_c^2) + brain_size_voxels +
        (1 | rat_code),
      data = datos_region,
      REML = FALSE
    ),
    error = function(e) NULL
  )
  
  
  # Guardar modelos
  modelos_lineales[[reg]] <- mod_linear
  modelos_cuadraticos[[reg]] <- mod_quad
  
  
  # ==========================
  # Comparación BIC
  # ==========================
  
  if (!is.null(mod_linear) && !is.null(mod_quad)) {
    
    bic_linear <- BIC(mod_linear)
    bic_quad <- BIC(mod_quad)
    
    resultados_BIC <- bind_rows(
      resultados_BIC,
      data.frame(
        region = reg,
        BIC_linear = bic_linear,
        BIC_quadratic = bic_quad,
        Delta_BIC = bic_quad - bic_linear,
        mejor_modelo = ifelse(
          bic_quad < bic_linear,
          "Quadratic",
          "Linear"
        )
      )
    )
  }
}

resultados_BIC
library(lme4)
library(lmerTest)
library(dplyr)

# ============================================================
# Regiones seleccionadas por BIC como modelo cuadrático
# ============================================================

regiones_cuadraticas <- resultados_BIC %>%
  filter(mejor_modelo == "Quadratic") %>%
  pull(region)


# ============================================================
# Ajustar modelo cuadrático y extraer estadística
# ============================================================

modelos_quad <- list()
resultados_quad <- data.frame()

for (reg in regiones_cuadraticas) {
  
  datos_region <- wm_data %>%
    filter(region == reg) %>%
    mutate(
      pnd_reg_c = pnd - mean(pnd, na.rm = TRUE)
    )
  
  # Modelo cuadrático
  mod_quad <- tryCatch(
    lmer(
      MD ~ pnd_reg_c * sex +
        I(pnd_reg_c^2) + brain_size_voxels +
        (1 | rat_code),
      data = datos_region,
      REML = TRUE
    ),
    error = function(e) NULL
  )
  
  # Guardar modelo
  modelos_quad[[reg]] <- mod_quad
  
  
  # ==========================================================
  # Extraer resultados
  # ==========================================================
  
  if (!is.null(mod_quad)) {
    
    coefs <- summary(mod_quad)$coefficients
    
    resultado <- data.frame(
      
      region = reg,
      
      # PND lineal
      estimate_pnd = coefs["pnd_reg_c", "Estimate"],
      SE_pnd = coefs["pnd_reg_c", "Std. Error"],
      t_pnd = coefs["pnd_reg_c", "t value"],
      p_pnd = coefs["pnd_reg_c", "Pr(>|t|)"],
      
      # PND cuadrático
      estimate_pnd2 = coefs["I(pnd_reg_c^2)", "Estimate"],
      SE_pnd2 = coefs["I(pnd_reg_c^2)", "Std. Error"],
      t_pnd2 = coefs["I(pnd_reg_c^2)", "t value"],
      p_pnd2 = coefs["I(pnd_reg_c^2)", "Pr(>|t|)"],
      
      # Sexo
      estimate_sex = coefs["sexMale", "Estimate"],
      SE_sex = coefs["sexMale", "Std. Error"],
      t_sex = coefs["sexMale", "t value"],
      p_sex = coefs["sexMale", "Pr(>|t|)"],
      
      # Interacción PND × sexo
      estimate_pnd_sex = coefs["pnd_reg_c:sexMale", "Estimate"],
      SE_pnd_sex = coefs["pnd_reg_c:sexMale", "Std. Error"],
      t_pnd_sex = coefs["pnd_reg_c:sexMale", "t value"],
      p_pnd_sex = coefs["pnd_reg_c:sexMale", "Pr(>|t|)"]
    )
    
    resultados_quad <- bind_rows(
      resultados_quad,
      resultado
    )
  }
}


# ============================================================
# FDR
# ============================================================

resultados_quad <- resultados_quad %>%
  mutate(
    
    # FDR para PND
    p_pnd_FDR = p.adjust(
      p_pnd,
      method = "fdr"
    ),
    
    # FDR para PND²
    p_pnd2_FDR = p.adjust(
      p_pnd2,
      method = "fdr"
    ),
    
    # FDR para sexo
    p_sex_FDR = p.adjust(
      p_sex,
      method = "fdr"
    ),
    
    # FDR para interacción PND × sexo
    p_pnd_sex_FDR = p.adjust(
      p_pnd_sex,
      method = "fdr"
    )
  )


# ============================================================
# Significancia
# ============================================================

resultados_quad <- resultados_quad %>%
  mutate(
    sig_pnd = p_pnd_FDR < 0.05,
    sig_pnd2 = p_pnd2_FDR < 0.05,
    sig_sex = p_sex_FDR < 0.05,
    sig_pnd_sex = p_pnd_sex_FDR < 0.05
  )


# ============================================================
# Ver resultados
# ============================================================

resultados_quad

library(writexl)

write_xlsx(
  resultados_quad,
  "resultados_modelos_cuadraticos_MD_GM.xlsx"
)

library(ggplot2)
library(patchwork)

plots_quad <- list()

for (reg in regiones_cuadraticas) {
  
  datos_region <- wm_data %>%
    filter(region == reg) %>%
    mutate(
      pnd_reg_c = pnd - mean(pnd, na.rm = TRUE)
    )
  
  modelo <- modelos_quad[[reg]]
  
  # ==========================================
  # Crear valores de PND para las predicciones
  # ==========================================
  
  newdata <- expand.grid(
    pnd = seq(
      min(datos_region$pnd, na.rm = TRUE),
      max(datos_region$pnd, na.rm = TRUE),
      length.out = 200
    ),
    sex = levels(datos_region$sex)
  )
  
  # Centrar PND exactamente igual que en el modelo
  media_pnd <- mean(datos_region$pnd, na.rm = TRUE)
  
  newdata$pnd_reg_c <- newdata$pnd - media_pnd
  
  # Predicción poblacional
  newdata$pred_MD <- predict(
    modelo,
    newdata = newdata,
    re.form = NA
  )
  
  
  # ==========================================
  # Plot
  # ==========================================
  
  p <- ggplot(
    datos_region,
    aes(
      x = pnd,
      y = MD,
      color = sex
    )
  ) +
    
    # Datos
    geom_point(
      alpha = 0.5,
      size = 2
    ) +
    
    # Curvas cuadráticas
    geom_line(
      data = newdata,
      aes(
        x = pnd,
        y = pred_MD,
        color = sex
      ),
      linewidth = 1.2
    ) +
    
    labs(
      title = paste("MD - Región:", reg),
      subtitle = "Modelo cuadrático",
      x = "PND",
      y = "MD",
      color = "Sex"
    ) +
    
    theme_classic() +
    
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  plots_quad[[reg]] <- p
}
library(patchwork)
# Mostrar todos los plots juntos
wrap_plots(plots_quad)
