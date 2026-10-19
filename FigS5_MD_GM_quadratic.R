#### QUADRATIC LINEAR MIXED MODEL - MD / GM ####

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

# Cargar resultados de la comparación BIC
resultados_regiones <- read_excel(
  "D:/Ageing/Ampliación_paper_reviews/WM_MRI/LMM_GM_QUADRATIC/MD/BIC_Lineal_Cuadratico.xlsx"
)

# Seleccionar únicamente las regiones cuyo modelo ganador es el cuadrático
regiones_cuadraticas <- resultados_regiones %>%
  filter(modelo == "Cuadrático") %>%
  pull(region)

# Ver las regiones
regiones_cuadraticas

# ============================================================
# 1. CARGAR DATOS
# ============================================================

ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/MD/CleanedMD_without_Outliers/gm_data_MD_without_outliers.xlsx"

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

gm_regiones <- regiones_cuadraticas

# # Comprobar regiones
# gm_regiones <- levels(gm_data$region)

# ============================================================
# 3. MODELO GLOBAL CON EFECTO ALEATORIO DE REGION
# ============================================================

##Interacción solo age*sex
mod_global_quad <- lmer(
  MD ~ pnd_c * sex + I(pnd_c^2) + 
    (1 | rat_code) +
    (1 | region),
  data = gm_data,
  REML = TRUE
)
summary(mod_global_quad)
anova(mod_global_quad)

##AÑADIENDO BRAIN SIZE
mod_global_quad <- lmer(
  MD ~ pnd_c * sex + I(pnd_c^2) + brain_size_voxels +
    (1 | rat_code) +
    (1 | region),
  data = gm_data,
  REML = TRUE
)

summary(mod_global_quad)
anova(mod_global_quad)


# ============================================================
# 4. MODELO CUADRÁTICO PARA CADA REGIÓN
# ============================================================

modelos <- list()
resultados <- list()

dir.create("LMM_GM_QUADRATIC_BIC_DELETE_2", showWarnings = FALSE)

for(reg in gm_regiones){
  
  # Filtramos y centramos la edad dentro de cada región
  datos_region <- gm_data %>%
    filter(region == reg) %>%
    mutate(
      pnd_reg_c = pnd - mean(pnd, na.rm = TRUE)
    )
  
  # ----------------------------------------------------------
  # Modelo cuadrático por región
  # ----------------------------------------------------------
  
  mod <- tryCatch(
    lmer(
      MD ~ pnd_reg_c * sex + I(pnd_reg_c^2) + brain_size_voxels +
        (1 | rat_code),
      data = datos_region,
      REML = TRUE
    ),
    error = function(e) NULL
  )
  
  # Guardar modelo
  modelos[[reg]] <- mod
  
  # ----------------------------------------------------------
  # Extraer resultados
  # ----------------------------------------------------------
  
  if(!is.null(mod)){
    
    an <- anova(mod, type = 3)
    
    resultados[[reg]] <- data.frame(
      region = reg,
      p_age_linear =
        an["pnd_reg_c", "Pr(>F)"],
      p_age_quadratic =
        an["I(pnd_reg_c^2)", "Pr(>F)"],
      p_sex =
        an["sex", "Pr(>F)"],
      p_age_linear_sex =
        an["pnd_reg_c:sex", "Pr(>F)"],
      p_age_quadratic_sex =
        an["I(pnd_reg_c^2):sex", "Pr(>F)"]
    )
    
    # --------------------------------------------------------
    # Guardar summary
    # --------------------------------------------------------
    
    capture.output(
      summary(mod),
      file = file.path(
        "LMM_GM_QUADRATIC_BIC_DELETE_2",
        paste0("Summary_", reg, ".txt")
      )
    )
  }
}

# ============================================================
# 5. UNIR RESULTADOS
# ============================================================

resultados <- bind_rows(resultados)

# ============================================================
# 6. CORRECCIÓN FDR
# ============================================================

resultados <- resultados %>%
  mutate(
    
    p_age_linear_FDR =
      p.adjust(p_age_linear, method = "BH"),
    
    p_age_quadratic_FDR =
      p.adjust(p_age_quadratic, method = "BH"),
    
    p_sex_FDR =
      p.adjust(p_sex, method = "BH"),
    
    p_age_linear_sex_FDR =
      p.adjust(p_age_linear_sex, method = "BH")
  )

# ============================================================
# 6B. REGIONES CON INTERACCIÓN AGE × SEX SIGNIFICATIVA
# ============================================================

regiones_significativas_age_sex <- resultados %>%
  filter(p_age_linear_sex_FDR < 0.05) %>%
  dplyr::select(
    region,
    p_age_linear_sex,
    p_age_linear_sex_FDR
  )

regiones_significativas_age_sex

# ============================================================
# 7. GUARDAR RESULTADOS
# ============================================================

write.csv(
  resultados,
  "LMM_results_GM_QUADRATIC_BIC.csv",
  row.names = FALSE
)

resultados


####GRAFICOS
library(ggplot2)
library(dplyr)

# Crear carpeta para las figuras
dir.create(
  "LMM_GM_QUADRATIC_BIC/Figures_delete",
  showWarnings = FALSE,
  recursive = TRUE
)

# Colores
colores_sex <- c(
  "Male" = "#77DD77",
  "Female" = "#9370DB"
)

# ----------------------------------------------------------
# Plot para cada región
# ----------------------------------------------------------

for(reg in gm_regiones){
  
  # Recuperar modelo de esa región
  mod <- modelos[[reg]]
  
  # Si el modelo no existe o no convergió, saltamos
  if(is.null(mod)) next
  
  # --------------------------------------------------------
  # Datos de la región
  # --------------------------------------------------------
  
  datos_region <- gm_data %>%
    filter(region == reg) %>%
    mutate(
      pnd_reg_c = pnd - mean(pnd, na.rm = TRUE)
    )
  
  # --------------------------------------------------------
  # Secuencia de edad para dibujar la curva
  # --------------------------------------------------------
  
  edad_seq <- seq(
    min(datos_region$pnd, na.rm = TRUE),
    max(datos_region$pnd, na.rm = TRUE),
    length.out = 200
  )
  
  # Media de PND utilizada para centrar el modelo
  media_pnd <- mean(
    datos_region$pnd,
    na.rm = TRUE
  )
  
  # --------------------------------------------------------
  # Dataset para predicciones
  # --------------------------------------------------------
  
  pred_data <- expand.grid(
    pnd = edad_seq,
    sex = levels(datos_region$sex)
  )
  
  pred_data <- pred_data %>%
    mutate(
      pnd_reg_c = pnd - media_pnd
    )
  
  # --------------------------------------------------------
  # Predicciones del modelo
  # No incluimos el efecto aleatorio de rat_code
  # --------------------------------------------------------
  
  pred_data$pred <- predict(
    mod,
    newdata = pred_data,
    re.form = NA
  )
  
  # --------------------------------------------------------
  # Gráfico
  # --------------------------------------------------------
  
  p <- ggplot(
    datos_region,
    aes(
      x = pnd,
      y = MD,
      colour = sex
    )
  ) +
    
    # Puntos individuales
    geom_point(
      alpha = 0.35,
      size = 1.5
    ) +
    
    # Trayectorias individuales de las ratas
    geom_line(
      aes(group = rat_code),
      alpha = 0.15
    ) +
    
    # Curvas cuadráticas estimadas
    geom_line(
      data = pred_data,
      aes(
        x = pnd,
        y = pred,
        colour = sex
      ),
      linewidth = 1.3
    ) +
    
    # Colores
    scale_colour_manual(
      values = colores_sex
    ) +
    
    # Tema
    theme_classic() +
    
    # Etiquetas
    labs(
      title = reg,
      subtitle = "Quadratic linear mixed model",
      x = "PND (age)",
      y = "FR",
      colour = "Sex"
    )
  
  # --------------------------------------------------------
  # Guardar PDF
  # --------------------------------------------------------
  
  ggsave(
    filename = paste0(
      "LMM_GM_QUADRATIC_BIC/Figures_delete/",
      reg,
      "_quadratic.pdf"
    ),
    plot = p,
    width = 7,
    height = 5
  )
  
  # --------------------------------------------------------
  # Guardar también PNG
  # --------------------------------------------------------
  
  ggsave(
    filename = paste0(
      "LMM_GM_QUADRATIC_BIC/Figures/",
      reg,
      "_quadratic.png"
    ),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  # Mostrar el gráfico en RStudio
  print(p)
}

 # ============================================================
# T-TEST MACHOS vs HEMBRAS EN CADA EDAD
# ============================================================

regiones_sig <- regiones_significativas_age_sex$region

datos_ttest <- gm_data %>%
  filter(region %in% regiones_sig)

resultados_ttest <- datos_ttest %>%
  group_by(region, pnd) %>%
  group_modify(~ {
    
    datos <- .x
    
    # --------------------------------------------------------
    # Eliminar NA de MD y sex
    # --------------------------------------------------------
    
    datos <- datos %>%
      filter(
        !is.na(MD),
        !is.na(sex)
      )
    
    # Número de observaciones por sexo
    n_female <- sum(datos$sex == "Female")
    n_male   <- sum(datos$sex == "Male")
    
    # --------------------------------------------------------
    # Si no hay al menos 2 observaciones en cada sexo,
    # no hacemos el t-test
    # --------------------------------------------------------
    
    if(n_female < 2 || n_male < 2){
      
      return(data.frame(
        n_female = n_female,
        n_male = n_male,
        mean_female = NA_real_,
        mean_male = NA_real_,
        difference = NA_real_,
        p_value = NA_real_
      ))
    }
    
    # --------------------------------------------------------
    # Welch t-test
    # --------------------------------------------------------
    
    test <- t.test(
      MD ~ sex,
      data = datos,
      var.equal = FALSE
    )
    
    # --------------------------------------------------------
    # Medias y diferencia Male - Female
    # --------------------------------------------------------
    
    mean_female <- mean(
      datos$MD[datos$sex == "Female"],
      na.rm = TRUE
    )
    
    mean_male <- mean(
      datos$MD[datos$sex == "Male"],
      na.rm = TRUE
    )
    
    data.frame(
      n_female = n_female,
      n_male = n_male,
      mean_female = mean_female,
      mean_male = mean_male,
      difference = mean_male - mean_female,
      p_value = test$p.value
    )
    
  }) %>%
  ungroup()


# ============================================================
# CORRECCIÓN FDR
# ============================================================

resultados_ttest <- resultados_ttest %>%
  group_by(pnd) %>%
  mutate(
    p_FDR = p.adjust(
      p_value,
      method = "BH"
    ),
    significativo_FDR = p_FDR < 0.05
  ) %>%
  ungroup()

# ============================================================
# TODOS LOS RESULTADOS
# ============================================================

resultados_ttest


# ============================================================
# SOLO COMPARACIONES SIGNIFICATIVAS
# ============================================================

resultados_ttest_sig <- resultados_ttest %>%
  filter(
    significativo_FDR
  )

resultados_ttest_sig


# ============================================================
# GUARDAR RESULTADOS
# ============================================================

write.csv(
  resultados_ttest,
  "Ttest_Male_vs_Female_by_age_GM.csv",
  row.names = FALSE
)

write.csv(
  resultados_ttest_sig,
  "Ttest_Male_vs_Female_by_age_GM_SIGNIFICANT.csv",
  row.names = FALSE
)

# ============================================================
# REPRESENTACIÓN: HEATMAP DE DIFERENCIAS SEXUALES
# ============================================================

library(ggplot2)
library(dplyr)

# Crear variable de significación
heatmap_data <- resultados_ttest %>%
  mutate(
    significativo = ifelse(p_FDR < 0.05, "Significativo", "No significativo"),
    
    # Dirección de la diferencia:
    # Male - Female
    direccion = case_when(
      p_FDR < 0.05 & difference > 0 ~ "Male > Female",
      p_FDR < 0.05 & difference < 0 ~ "Female > Male",
      TRUE ~ "No significativo"
    )
  )

# Ordenar las edades
heatmap_data$pnd <- factor(
  heatmap_data$pnd,
  levels = sort(unique(heatmap_data$pnd))
)

# Ordenar regiones
heatmap_data$region <- factor(
  heatmap_data$region,
  levels = rev(unique(heatmap_data$region))
)

# ------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------

p_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = pnd,
    y = region,
    fill = direccion
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  
  scale_fill_manual(
    values = c(
      "Male > Female" = "#77DD77",
      "Female > Male" = "#9370DB",
      "No significativo" = "#EAEAEA"
    )
  ) +
  
  theme_classic() +
  
  labs(
    x = "PND",
    y = "Region",
    fill = "Sex difference"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )

print(p_heatmap)

ggsave(
  "TTEST_sex_by_age_heatmap.pdf",
  p_heatmap,
  width = 9,
  height = 6
)

