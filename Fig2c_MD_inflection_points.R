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
# 2. LISTAS PARA GUARDAR MODELOS Y VÉRTICES
# ============================================================

modelos <- list()
resultados <- list()

dir.create(
  "LMM_GM_QUADRATIC_BIC_figure2_1",
  showWarnings = FALSE
)


# ============================================================
# 3. REGIONES
# ============================================================

gm_regiones <- levels(gm_data$region)


# ============================================================
# 4. LOOP POR REGIÓN
# ============================================================

for(reg in gm_regiones){
  
  # ----------------------------------------------------------
  # Datos de la región
  # ----------------------------------------------------------
  
  datos_region <- gm_data %>%
    filter(region == reg) %>%
    mutate(
      pnd_reg_c = pnd - mean(pnd, na.rm = TRUE)
    )
  
  # Media de pnd de la región
  mean_pnd <- mean(
    datos_region$pnd,
    na.rm = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Rango de pnd observado
  # ----------------------------------------------------------
  
  min_pnd <- min(
    datos_region$pnd,
    na.rm = TRUE
  )
  
  max_pnd <- max(
    datos_region$pnd,
    na.rm = TRUE
  )
  
  
  # ==========================================================
  # MODELO CUADRÁTICO
  # ==========================================================
  
  mod <- tryCatch(
    
    lmer(
      MD ~ pnd_reg_c * sex +
        I(pnd_reg_c^2) +
        (1 | rat_code),
      data = datos_region,
      REML = TRUE
    ),
    
    error = function(e) NULL
  )
  
  
  # Guardar modelo
  modelos[[reg]] <- mod
  
  
  # ==========================================================
  # SI EL MODELO SE AJUSTÓ CORRECTAMENTE
  # ==========================================================
  
  if(!is.null(mod)){
    
    
    # --------------------------------------------------------
    # Coeficientes
    # --------------------------------------------------------
    
    coefs <- fixef(mod)
    
    beta_linear <- 
      unname(coefs["pnd_reg_c"])
    
    beta_quadratic <- 
      unname(coefs["I(pnd_reg_c^2)"])
    
    beta_age_sex <- 
      unname(coefs["pnd_reg_c:sexMale"])
    
    
    # ========================================================
    # VÉRTICE FEMALE
    # ========================================================
    
    vertex_female_c <- 
      -beta_linear /
      (2 * beta_quadratic)
    
    vertex_female <- 
      vertex_female_c +
      mean_pnd
    
    
    # ========================================================
    # VÉRTICE MALE
    # ========================================================
    
    beta_linear_male <- 
      beta_linear +
      beta_age_sex
    
    vertex_male_c <- 
      -beta_linear_male /
      (2 * beta_quadratic)
    
    vertex_male <- 
      vertex_male_c +
      mean_pnd
    
    
    # ========================================================
    # MÍNIMO O MÁXIMO
    # ========================================================
    
    vertex_type <- ifelse(
      beta_quadratic > 0,
      "minimum",
      "maximum"
    )
    
    
    # ========================================================
    # ¿ESTÁ DENTRO DEL RANGO OBSERVADO?
    # ========================================================
    
    vertex_female_in_range <-
      vertex_female >= min_pnd &
      vertex_female <= max_pnd
    
    vertex_male_in_range <-
      vertex_male >= min_pnd &
      vertex_male <= max_pnd
    
    
    # ========================================================
    # GUARDAR RESULTADOS
    # ========================================================
    
    resultados[[reg]] <- data.frame(
      
      region = reg,
      
      # Coeficientes
      beta_linear = beta_linear,
      beta_quadratic = beta_quadratic,
      beta_age_sex = beta_age_sex,
      
      # Vértices en escala original
      vertex_female = vertex_female,
      vertex_male = vertex_male,
      
      # Vértices en escala centrada
      vertex_female_centered =
        vertex_female_c,
      
      vertex_male_centered =
        vertex_male_c,
      
      # Tipo de vértice
      vertex_type = vertex_type,
      
      # ¿Está dentro del rango observado?
      vertex_female_in_range =
        vertex_female_in_range,
      
      vertex_male_in_range =
        vertex_male_in_range,
      
      # Información sobre pnd
      min_pnd = min_pnd,
      max_pnd = max_pnd,
      mean_pnd = mean_pnd
    )
    
    
    # ========================================================
    # GUARDAR SUMMARY DEL MODELO
    # ========================================================
    
    capture.output(
      summary(mod),
      file = file.path(
        "LMM_GM_QUADRATIC_BIC_figure2_1",
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
# 6. GUARDAR RESULTADOS
# ============================================================

write.csv(
  resultados,
  "LMM_results_GM_QUADRATIC_vertices.csv",
  row.names = FALSE
)

####GRAFICOS
library(ggplot2)
library(dplyr)

# Crear carpeta para las figuras
dir.create(
  "LMM_GM_QUADRATIC_BIC_figure2_!/Figures",
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
      y = "MD",
      colour = "Sex"
    )
  
  # --------------------------------------------------------
  # Guardar PDF
  # --------------------------------------------------------
  
  ggsave(
    filename = paste0(
      "LMM_GM_QUADRATIC_BIC_figure2_1/Figures/",
      reg,
      "_quadratic.pdf"
    ),
    plot = p,
    width = 7,
    height = 5
  )
}
