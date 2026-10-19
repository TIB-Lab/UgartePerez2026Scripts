rm(list=ls())
library(lme4)
library(lmerTest)
library(emmeans)

# Especificar la ruta del archivo xlsx
#FA:
#ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/FA/NoCleaned_FA_without_outliers/wm_data_FA_without_outliers.xlsx"
#FR:
ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/FR/WM/wm_data_FR_without_outliers.xlsx"
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

modelo_FA_lineal <- lmer(
  FR ~ pnd * sex  + brain_size_voxels + 
    (1 | rat_code) + 
    (1 | region),
  data = wm_data,
  REML = FALSE
)
modelo_FA_cuadratico <- lmer(
    FA ~ pnd * sex  + brain_size_voxels +
      I(pnd^2) +
      (1 | rat_code)+
      (1 | region),
    data = wm_data,
    REML= FALSE 
  )

# Comparar BIC BIC(modelo_FA_lineal, modelo_FA_cuadratico)
BIC(modelo_FA_lineal, modelo_FA_cuadratico)

summary(modelo_FA_cuadratico)
anova(modelo_FA_cuadratico)

# 
# wm_regiones <- setdiff(
#   wm_regiones,
#   "59_WM_fas_retrof"
# )

modelo <- lmer(
  FA ~
    pnd * sex + I(pnd^2) + brain_size_voxels +
    (1 | rat_code) +
    (1 | region),
  data = wm_data,
  REML = TRUE
)

summary(modelo)
anova(modelo)


library(lme4)
library(lmerTest)
library(dplyr)

# ==========================================
# 1. Ajustar un modelo para cada región
# ==========================================

modelos_FA <- list()

for (reg in wm_regiones) {
  
  datos_region <- wm_data %>%
    filter(region == reg)
  
  modelos_FA[[reg]] <- lmer(
    
    FR ~ pnd * sex + I(pnd^2) + brain_size_voxels +
      (1 | rat_code),
    data = datos_region
  )
}


# ==========================================
# 2. Extraer sex y pnd:sex
# ==========================================

resultados <- data.frame()

for (reg in names(modelos_FA)) {
  
  modelo <- modelos_FA[[reg]]
  anova_modelo <- anova(modelo)
  
  resultado <- data.frame(
    region = reg,
    
    # Efecto principal de sex
    p_sex = anova_modelo["sex", "Pr(>F)"],
    
    # Interacción pnd × sex
    p_pnd_sex = anova_modelo["pnd:sex", "Pr(>F)"]
  )
  
  resultados <- bind_rows(
    resultados,
    resultado
  )
}


# ==========================================
# 3. Corrección FDR
# ==========================================

resultados <- resultados %>%
  mutate(
    
    # FDR para sex entre todas las regiones
    p_sex_FDR = p.adjust(
      p_sex,
      method = "fdr"
    ),
    
    # FDR para pnd × sex entre todas las regiones
    p_pnd_sex_FDR = p.adjust(
      p_pnd_sex,
      method = "fdr"
    )
  )


# ==========================================
# 4. Significancia
# ==========================================

resultados <- resultados %>%
  mutate(
    sex_significativo = p_sex_FDR < 0.05,
    interaccion_significativa = p_pnd_sex_FDR < 0.05
  )


# Ver resultados
resultados

write.csv(
  resultados,
  "LMM_results_FA_WM_QUADRATIC_BIC_brainsize.csv",
  row.names = FALSE
)


###PLOT 
library(ggplot2)

# Crear un plot para cada región
plots_FA <- list()

for (reg in wm_regiones) {
  
  datos_region <- wm_data %>%
    filter(region == reg)
  
  modelo <- modelos_FA[[reg]]
  
  # Predicciones del modelo
  datos_region$pred_FA <- predict(modelo, re.form = NA)
  
  # Plot
  p <- ggplot(datos_region, aes(x = pnd, y = MD, color = sex)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_line(
      aes(y = pred_MD),
      linewidth = 1
    ) +
    labs(
      title = paste("FA - Región:", reg),
      x = "PND",
      y = "FA",
      color = "Sex"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  plots_FA[[reg]] <- p
}
library(patchwork)
# Mostrar todos los plots juntos
wrap_plots(plots_FA)
