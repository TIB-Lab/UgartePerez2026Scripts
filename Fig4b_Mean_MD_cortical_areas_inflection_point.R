# ============================================================
# LMM CUADRÁTICO:
# MEDIA DE MD EN REGIONES CON AGE × SEX SIGNIFICATIVO TRAS FDR
#
# OBJETIVOS:
# 1. Seleccionar las regiones significativas tras FDR.
# 2. Calcular la media de MD por rata × PND × sexo.
# 3. Ajustar un único LMM cuadrático.
# 4. Evaluar PND, PND², sexo y PND × sexo.
# 5. Obtener las trayectorias predichas para ambos sexos.
# 6. Calcular los mínimos de las trayectorias.
# 7. Calcular IC95% de la media ajustada por efectos fijos.
# 8. Representar y exportar todos los resultados.
#
# MODELO:
# MD_mean_regions ~ pnd_c * sex + I(pnd_c^2) + (1 | rat_code)
#
# ============================================================


# ============================================================
# 0. LIBRERÍAS
# ============================================================

library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(openxlsx)


# ============================================================
# 1. CONFIGURACIÓN
# ============================================================

id_rata <- "rat_code"

carpeta_salida <- file.path(
  "LMM_GM_QUADRATIC_BIC",
  "Figures"
)

dir.create(
  carpeta_salida,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. COMPROBAR QUE EXISTEN LOS OBJETOS Y LAS VARIABLES
# ============================================================

if (!exists("gm_data")) {
  stop("No existe el objeto gm_data.")
}

if (!exists("regiones_significativas_age_sex")) {
  stop(
    "No existe el objeto regiones_significativas_age_sex."
  )
}

variables_necesarias <- c(
  id_rata,
  "region",
  "pnd",
  "sex",
  "MD"
)

faltan <- setdiff(
  variables_necesarias,
  names(gm_data)
)

if (length(faltan) > 0) {
  stop(
    "Faltan columnas en gm_data: ",
    paste(faltan, collapse = ", ")
  )
}

if (!"region" %in% names(regiones_significativas_age_sex)) {
  stop(
    "El objeto regiones_significativas_age_sex ",
    "debe contener una columna llamada region."
  )
}


# ============================================================
# 3. SELECCIONAR REGIONES SIGNIFICATIVAS TRAS FDR
# ============================================================

regiones_sig <- unique(
  as.character(
    regiones_significativas_age_sex$region
  )
)

regiones_sig <- regiones_sig[
  !is.na(regiones_sig)
]

if (length(regiones_sig) == 0) {
  stop(
    "No se han encontrado regiones significativas."
  )
}

regiones_disponibles <- intersect(
  regiones_sig,
  unique(as.character(gm_data$region))
)

regiones_no_disponibles <- setdiff(
  regiones_sig,
  regiones_disponibles
)

if (length(regiones_no_disponibles) > 0) {
  warning(
    "Estas regiones no aparecen en gm_data: ",
    paste(regiones_no_disponibles, collapse = ", ")
  )
}

regiones_sig <- regiones_disponibles

if (length(regiones_sig) == 0) {
  stop(
    "Ninguna región seleccionada aparece en gm_data."
  )
}

cat("\n============================================\n")
cat("REGIONES SIGNIFICATIVAS TRAS FDR\n")
cat("============================================\n")

print(regiones_sig)

cat(
  "\nNúmero de regiones seleccionadas:",
  length(regiones_sig),
  "\n"
)


# ============================================================
# 4. CREAR MEDIA DE MD POR RATA × PND × SEXO
# ============================================================

# Cada observación representa una rata en un PND.
# MD_mean_regions es la media de las regiones seleccionadas
# que tienen MD disponible para esa observación.

datos_individuales <- gm_data %>%
  filter(
    as.character(region) %in% regiones_sig
  ) %>%
  group_by(
    .data[[id_rata]],
    pnd,
    sex
  ) %>%
  summarise(
    MD_mean_regions = if (all(is.na(MD))) {
      NA_real_
    } else {
      mean(MD, na.rm = TRUE)
    },
    
    n_regions = sum(!is.na(MD)),
    
    .groups = "drop"
  ) %>%
  filter(
    n_regions > 0,
    !is.na(pnd),
    !is.na(sex),
    !is.na(.data[[id_rata]])
  )

if (nrow(datos_individuales) == 0) {
  stop(
    "No hay observaciones válidas después de calcular ",
    "la media de MD."
  )
}

# Asegurar que el sexo es un factor
datos_individuales$sex <- droplevels(
  factor(datos_individuales$sex)
)

niveles_sex <- levels(datos_individuales$sex)

if (length(niveles_sex) != 2) {
  stop(
    "El análisis requiere exactamente dos niveles de sex. ",
    "Niveles encontrados: ",
    paste(niveles_sex, collapse = ", ")
  )
}

cat("\n============================================\n")
cat("OBSERVACIONES POR SEXO\n")
cat("============================================\n")

print(
  table(datos_individuales$sex)
)

cat(
  "\nNúmero de observaciones rata × PND:",
  nrow(datos_individuales),
  "\n"
)

cat(
  "\nNúmero de ratas:",
  dplyr::n_distinct(datos_individuales[[id_rata]]),
  "\n"
)

cat(
  "\nNúmero de regiones disponibles:",
  length(regiones_sig),
  "\n"
)

# Revisar el número de regiones que contribuyen a cada media
print(
  summary(datos_individuales$n_regions)
)


# ============================================================
# 5. CENTRAR PND
# ============================================================

media_pnd <- mean(
  datos_individuales$pnd,
  na.rm = TRUE
)

datos_media_modelo <- datos_individuales %>%
  mutate(
    pnd_c = pnd - media_pnd
  )

cat(
  "\nMedia global de PND utilizada para centrar:",
  media_pnd,
  "\n"
)


# ============================================================
# 6. CONFIGURAR CONTRASTES PARA ANOVA TIPO III
# ============================================================

# Contrastes de suma para el factor sex.
# Esto permite interpretar adecuadamente el ANOVA tipo III
# con la codificación elegida.

contrasts(datos_media_modelo$sex) <- contr.sum(2)

# Comprobar que el identificador de rata existe
if (!id_rata %in% names(datos_media_modelo)) {
  stop(
    "No se encuentra el identificador de rata."
  )
}


# ============================================================
# 7. AJUSTAR UN ÚNICO LMM CUADRÁTICO
# ============================================================

# Modelo:
#
# MD_mean_regions ~ pnd_c * sex +
#                   I(pnd_c^2) +
#                   (1 | rat_code)
#
# PND² es común a ambos sexos.
# PND × sexo permite diferencias en el componente lineal.
# Se incluye un intercepto aleatorio por rata.

formula_modelo <- as.formula(
  paste0(
    "MD_mean_regions ~ ",
    "pnd_c * sex + ",
    "I(pnd_c^2) + ",
    "(1 | ", id_rata, ")"
  )
)

modelo_media <- lmer(
  formula_modelo,
  data = datos_media_modelo,
  REML = TRUE
)

cat("\n============================================\n")
cat("RESUMEN DEL LMM CUADRÁTICO\n")
cat("============================================\n")

print(
  summary(modelo_media)
)


# ============================================================
# 8. ANOVA TIPO III
# ============================================================

anova_media <- anova(
  modelo_media,
  type = 3
)

cat("\n============================================\n")
cat("ANOVA TIPO III\n")
cat("============================================\n")

print(anova_media)

# Interpretación de los términos:
#
# pnd_c:
#   Efecto lineal de PND bajo la codificación utilizada.
#
# I(pnd_c^2):
#   Efecto cuadrático común a ambos sexos.
#
# sex:
#   Diferencia entre sexos en PND centrado = 0.
#   Con contrastes de suma, representa el contraste
#   de sexo en la edad media de la muestra.
#
# pnd_c:sex:
#   Diferencia entre sexos en el componente lineal.
#
# Nota:
# El contraste lineal de PND depende de la codificación
# de sex. No es un efecto lineal promedio independiente
# de la parametrización.


# ============================================================
# 9. GUARDAR COEFICIENTES Y ANOVA
# ============================================================

beta <- lme4::fixef(modelo_media)

coeficientes <- as.data.frame(
  coef(summary(modelo_media))
)

coeficientes$Termino <- rownames(coeficientes)

rownames(coeficientes) <- NULL

coeficientes <- coeficientes %>%
  dplyr::select(
    Termino,
    everything()
  )

anova_tabla <- as.data.frame(anova_media)

anova_tabla$Termino <- rownames(anova_tabla)

rownames(anova_tabla) <- NULL

anova_tabla <- anova_tabla %>%
  dplyr::select(
    Termino,
    everything()
  )


# ============================================================
# 10. GENERAR SECUENCIA DE PND PARA LAS PREDICCIONES
# ============================================================

edad_seq <- seq(
  from = min(
    datos_media_modelo$pnd,
    na.rm = TRUE
  ),
  to = max(
    datos_media_modelo$pnd,
    na.rm = TRUE
  ),
  length.out = 200
)


# ============================================================
# 11. CREAR DATASET DE PREDICCIONES
# ============================================================

pred_data <- expand.grid(
  pnd = edad_seq,
  sex = niveles_sex,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

pred_data$sex <- factor(
  pred_data$sex,
  levels = niveles_sex
)

pred_data$pnd_c <- pred_data$pnd - media_pnd

# Añadir un identificador existente por si la función
# de predicción requiere la variable del efecto aleatorio.
# re.form = NA excluirá el efecto aleatorio de la predicción.

pred_data[[id_rata]] <- datos_media_modelo[[id_rata]][1]


# ============================================================
# 12. PREDICCIONES DEL LMM
# ============================================================

# Predicciones marginales basadas exclusivamente
# en los efectos fijos del modelo.

pred_data$pred <- predict(
  modelo_media,
  newdata = pred_data,
  re.form = NA
)


# ============================================================
# 13. IC95% DEL FIT DEL LMM
# ============================================================

# Se construye la matriz de diseño de efectos fijos
# con la misma fórmula y los mismos contrastes del modelo.

formula_fija <- lme4::nobars(
  formula(modelo_media)
)

terminos_fijos <- delete.response(
  terms(formula_fija)
)

# Recuperar los contrastes utilizados por el modelo
X_modelo <- model.matrix(modelo_media)

contrastes_modelo <- attr(
  X_modelo,
  "contrasts"
)

# Crear la matriz de diseño de las predicciones
X <- model.matrix(
  terminos_fijos,
  data = pred_data,
  contrasts.arg = contrastes_modelo
)

beta <- lme4::fixef(modelo_media)

V <- as.matrix(
  vcov(modelo_media)
)

# Comprobar correspondencia de coeficientes y columnas
if (!all(names(beta) %in% colnames(X))) {
  
  cat("\nNombres de los coeficientes:\n")
  print(names(beta))
  
  cat("\nColumnas de la matriz X:\n")
  print(colnames(X))
  
  stop(
    "Las columnas de la matriz X no coinciden ",
    "con los coeficientes del modelo."
  )
}

# Ordenar X y V según los coeficientes del modelo
X <- X[
  ,
  names(beta),
  drop = FALSE
]

V <- V[
  names(beta),
  names(beta),
  drop = FALSE
]

# Varianza de cada media predicha:
# Var(y_hat) = X %*% V %*% t(X)
#
# El cálculo siguiente obtiene la diagonal de esa matriz.

pred_data$se_pred <- sqrt(
  pmax(
    0,
    rowSums(
      (X %*% V) * X
    )
  )
)

# IC95% aproximado de la media ajustada
pred_data <- pred_data %>%
  mutate(
    CI95_fit_low =
      pred - 1.96 * se_pred,
    
    CI95_fit_high =
      pred + 1.96 * se_pred
  )

cat("\n============================================\n")
cat("COMPROBACIÓN DE LAS PREDICCIONES E IC95%\n")
cat("============================================\n")

print(
  head(
    pred_data[
      ,
      c(
        "pnd",
        "sex",
        "pred",
        "se_pred",
        "CI95_fit_low",
        "CI95_fit_high"
      )
    ]
  )
)


# ============================================================
# 14. CALCULAR EL MÍNIMO DE CADA TRAYECTORIA
# ============================================================

# Para una trayectoria:
#
# y = a + b*x + c*x²
#
# el vértice se encuentra en:
#
# x_min = -b / (2*c)
#
# Se obtienen los coeficientes de cada parábola a partir
# de las predicciones del MISMO modelo en x = -1, 0 y 1.
#
# Así evitamos depender de nombres de interacción que
# cambian según la codificación del factor sex.

calcular_minimo <- function(s) {
  
  nd <- data.frame(
    pnd_c = c(-1, 0, 1),
    sex = factor(
      s,
      levels = niveles_sex
    )
  )
  
  # Añadir el identificador de rata
  nd[[id_rata]] <- datos_media_modelo[[id_rata]][1]
  
  # Predicciones de efectos fijos en los tres puntos
  y <- as.numeric(
    predict(
      modelo_media,
      newdata = nd,
      re.form = NA
    )
  )
  
  # Coeficiente lineal específico de este sexo
  b_linear_s <- (y[3] - y[1]) / 2
  
  # Coeficiente cuadrático específico de esta trayectoria
  b_quad_s <- (y[3] + y[1] - 2 * y[2]) / 2
  
  # Comprobar si existe un mínimo parabólico
  if (
    !is.finite(b_quad_s) ||
    b_quad_s <= 0
  ) {
    
    return(
      data.frame(
        sex = s,
        PND_min = NA_real_,
        MD_min = NA_real_,
        dentro_rango = FALSE
      )
    )
  }
  
  # Mínimo en PND centrado
  pnd_min_c <- -b_linear_s / (2 * b_quad_s)
  
  # Volver a la escala original de PND
  pnd_min <- pnd_min_c + media_pnd
  
  # Obtener MD en el mínimo
  nd_min <- data.frame(
    pnd_c = pnd_min_c,
    sex = factor(
      s,
      levels = niveles_sex
    )
  )
  
  nd_min[[id_rata]] <- datos_media_modelo[[id_rata]][1]
  
  md_min <- as.numeric(
    predict(
      modelo_media,
      newdata = nd_min,
      re.form = NA
    )
  )
  
  data.frame(
    sex = s,
    PND_min = pnd_min,
    MD_min = md_min,
    dentro_rango =
      pnd_min >= min(datos_media_modelo$pnd, na.rm = TRUE) &&
      pnd_min <= max(datos_media_modelo$pnd, na.rm = TRUE)
  )
}

# Calcular mínimos de ambos sexos
minimos_media <- do.call(
  rbind,
  lapply(
    niveles_sex,
    calcular_minimo
  )
)

rownames(minimos_media) <- NULL

cat("\n============================================\n")
cat("MÍNIMOS DE LAS TRAYECTORIAS\n")
cat("============================================\n")

print(minimos_media)


# ============================================================
# 15. DEFINIR COLORES POR SEXO
# ============================================================

# Colores por defecto para Male y Female.
# Si los niveles tienen otros nombres, se asignan
# los mismos dos colores en el orden de los niveles.

colores_predefinidos <- c(
  "Male" = "#77DD77",
  "Female" = "#9370DB"
)

if (all(niveles_sex %in% names(colores_predefinidos))) {
  
  colores_sex <- colores_predefinidos[niveles_sex]
  
} else {
  
  colores_sex <- setNames(
    c("#77DD77", "#9370DB"),
    niveles_sex
  )
}


# ============================================================
# 16. REPRESENTACIÓN DEL MODELO CUADRÁTICO
# ============================================================

p_media <- ggplot() +
  
  # ----------------------------------------------------------
# IC95% DEL FIT
# ----------------------------------------------------------

geom_ribbon(
  data = pred_data,
  aes(
    x = pnd,
    ymin = CI95_fit_low,
    ymax = CI95_fit_high,
    fill = sex,
    group = sex
  ),
  alpha = 0.22,
  colour = NA
) +
  
  # ----------------------------------------------------------
# PUNTOS INDIVIDUALES
# Cada punto = rata × PND
# ----------------------------------------------------------

geom_point(
  data = datos_individuales,
  aes(
    x = pnd,
    y = MD_mean_regions,
    colour = sex
  ),
  alpha = 0.40,
  size = 1.8
) +
  
  # ----------------------------------------------------------
# TRAYECTORIAS CUADRÁTICAS DEL MISMO LMM
# ----------------------------------------------------------

geom_line(
  data = pred_data,
  aes(
    x = pnd,
    y = pred,
    colour = sex,
    group = sex
  ),
  linewidth = 1.3
) +
  
  # ----------------------------------------------------------
# MÍNIMOS DENTRO DEL RANGO OBSERVADO
# ----------------------------------------------------------

geom_point(
  data = minimos_media %>%
    filter(dentro_rango),
  aes(
    x = PND_min,
    y = MD_min,
    colour = sex
  ),
  size = 3.5
) +
  
  # ----------------------------------------------------------
# LÍNEAS VERTICALES EN LOS MÍNIMOS
# ----------------------------------------------------------

geom_vline(
  data = minimos_media %>%
    filter(dentro_rango),
  aes(
    xintercept = PND_min,
    colour = sex
  ),
  linetype = "dashed",
  linewidth = 0.6,
  alpha = 0.7
) +
  
  # ----------------------------------------------------------
# ETIQUETAS DE LOS MÍNIMOS
# ----------------------------------------------------------

geom_text(
  data = minimos_media %>%
    filter(dentro_rango),
  aes(
    x = PND_min,
    y = MD_min,
    label = paste0(
      "PND = ",
      round(PND_min, 1),
      "\nMD = ",
      round(MD_min, 4)
    ),
    colour = sex
  ),
  vjust = -1,
  size = 3.3,
  fontface = "bold",
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# ESCALAS
# ----------------------------------------------------------

scale_colour_manual(
  values = colores_sex,
  breaks = niveles_sex
) +
  
  scale_fill_manual(
    values = colores_sex,
    breaks = niveles_sex
  ) +
  
  # ----------------------------------------------------------
# TEMA Y ETIQUETAS
# ----------------------------------------------------------

theme_classic(
  base_size = 12
) +
  
  labs(
    title = "Quadratic mixed-model trajectories",
    
    subtitle = paste0(
      "Mean MD across regions with significant Age × Sex ",
      "interaction after FDR (n = ",
      length(regiones_sig),
      " regions)"
    ),
    
    x = "Postnatal day (PND)",
    
    y = "Mean MD across significant regions",
    
    colour = "Sex",
    
    fill = "Sex"
  )

# Mostrar gráfico
print(p_media)


# ============================================================
# 17. EXPORTAR FIGURA PDF
# ============================================================

ggsave(
  filename = file.path(
    carpeta_salida,
    "Mean_AgeSex_significant_regions_LMM.pdf"
  ),
  plot = p_media,
  width = 8,
  height = 5.5
)


# ============================================================
# 18. EXPORTAR FIGURA PNG
# ============================================================

ggsave(
  filename = file.path(
    carpeta_salida,
    "Mean_AgeSex_significant_regions_LMM.png"
  ),
  plot = p_media,
  width = 8,
  height = 5.5,
  dpi = 300
)


# ============================================================
# 19. EXPORTAR RESULTADOS DEL MODELO A EXCEL
# ============================================================

resultados_modelo <- list(
  
  Coeficientes = coeficientes,
  
  ANOVA_tipo_III = anova_tabla,
  
  Minimos = minimos_media,
  
  Regiones_significativas = data.frame(
    region = regiones_sig
  )
)

openxlsx::write.xlsx(
  resultados_modelo,
  file = file.path(
    carpeta_salida,
    "Mean_AgeSex_LMM_results.xlsx"
  ),
  overwrite = TRUE
)


# ============================================================
# 20. EXPORTAR MÍNIMOS
# ============================================================

write.csv(
  minimos_media,
  file = file.path(
    carpeta_salida,
    "Minimum_PND_mean_regions_LMM.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. EXPORTAR DATOS INDIVIDUALES
# ============================================================

write.csv(
  datos_individuales,
  file = file.path(
    carpeta_salida,
    "Individual_mean_regions_AgeSex.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. EXPORTAR PREDICCIONES E IC95%
# ============================================================

write.csv(
  pred_data,
  file = file.path(
    carpeta_salida,
    "Predicted_mean_regions_AgeSex.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. EXPORTAR RESUMEN COMPLETO DEL MODELO
# ============================================================

capture.output(
  summary(modelo_media),
  file = file.path(
    carpeta_salida,
    "Summary_Mean_AgeSex_LMM.txt"
  )
)


# ============================================================
# 24. EXPORTAR ANOVA
# ============================================================

capture.output(
  anova_media,
  file = file.path(
    carpeta_salida,
    "ANOVA_Type_III_Mean_AgeSex_LMM.txt"
  )
)


# ============================================================
# 25. RESUMEN FINAL
# ============================================================

cat("\n============================================\n")
cat("ANÁLISIS COMPLETADO\n")
cat("============================================\n")

cat(
  "\nModelo ajustado:\n"
)

print(
  formula(modelo_media)
)

cat(
  "\nNúmero de regiones:",
  length(regiones_sig),
  "\n"
)

cat(
  "Número de observaciones:",
  nrow(datos_media_modelo),
  "\n"
)

cat(
  "Número de ratas:",
  dplyr::n_distinct(datos_media_modelo[[id_rata]]),
  "\n"
)

cat(
  "\nMínimos estimados:\n"
)

print(
  minimos_media
)

cat(
  "\nArchivos guardados en:\n",
  normalizePath(
    carpeta_salida,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)

# ============================================================
# REPRESENTACIÓN DEL LMM CUADRÁTICO
# ============================================================

library(ggplot2)
library(dplyr)
library(lme4)

# ------------------------------------------------------------
# 1. GENERAR PREDICCIONES PARA AMBOS SEXOS
# ------------------------------------------------------------

niveles_sex <- levels(datos_media_modelo$sex)

edad_seq <- seq(
  min(datos_media_modelo$pnd, na.rm = TRUE),
  max(datos_media_modelo$pnd, na.rm = TRUE),
  length.out = 200
)

pred_data <- expand.grid(
  pnd = edad_seq,
  sex = niveles_sex,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

pred_data$sex <- factor(
  pred_data$sex,
  levels = niveles_sex
)

pred_data$pnd_c <- pred_data$pnd - media_pnd

# Identificador válido para predict()
pred_data[[id_rata]] <- datos_media_modelo[[id_rata]][1]

# ------------------------------------------------------------
# 2. PREDICCIONES DEL MODELO
# ------------------------------------------------------------

pred_data$pred <- predict(
  modelo_media,
  newdata = pred_data,
  re.form = NA
)

# ------------------------------------------------------------
# 3. IC95% DEL FIT A PARTIR DE LOS EFECTOS FIJOS
# ------------------------------------------------------------

formula_fija <- lme4::nobars(
  formula(modelo_media)
)

terminos_fijos <- delete.response(
  terms(formula_fija)
)

# Recuperar contrastes utilizados en el modelo
X_modelo <- model.matrix(modelo_media)

contrastes_modelo <- attr(
  X_modelo,
  "contrasts"
)

X <- model.matrix(
  terminos_fijos,
  data = pred_data,
  contrasts.arg = contrastes_modelo
)

beta <- lme4::fixef(modelo_media)
V <- as.matrix(vcov(modelo_media))

# Comprobar correspondencia entre columnas y coeficientes
if (!all(names(beta) %in% colnames(X))) {
  stop(
    "Las columnas de X no coinciden con los coeficientes. ",
    "Revisa names(fixef(modelo_media)) y colnames(X)."
  )
}

X <- X[, names(beta), drop = FALSE]
V <- V[names(beta), names(beta), drop = FALSE]

pred_data$se_pred <- sqrt(
  pmax(
    0,
    rowSums((X %*% V) * X)
  )
)

pred_data <- pred_data %>%
  mutate(
    CI95_fit_low = pred - 1.96 * se_pred,
    CI95_fit_high = pred + 1.96 * se_pred
  )

# ------------------------------------------------------------
# 4. COLORES
# ------------------------------------------------------------

colores_base <- c(
  "Male" = "#77DD77",
  "Female" = "#9370DB"
)

if (all(niveles_sex %in% names(colores_base))) {
  colores_sex <- colores_base[niveles_sex]
} else {
  colores_sex <- setNames(
    c("#77DD77", "#9370DB"),
    niveles_sex
  )
}

# ------------------------------------------------------------
# 5. PREPARAR LOS MÍNIMOS
# ------------------------------------------------------------

minimos_plot <- minimos_media %>%
  filter(
    dentro_rango,
    !is.na(PND_min),
    !is.na(MD_min)
  )

# ------------------------------------------------------------
# 6. CONSTRUIR EL GRÁFICO
# ------------------------------------------------------------

p_media <- ggplot() +
  
  # IC95% de la media ajustada
  geom_ribbon(
    data = pred_data,
    aes(
      x = pnd,
      ymin = CI95_fit_low,
      ymax = CI95_fit_high,
      fill = sex,
      group = sex
    ),
    alpha = 0.22,
    colour = NA
  ) +
  
  # Observaciones individuales:
  # cada punto representa una rata × PND
  geom_point(
    data = datos_individuales,
    aes(
      x = pnd,
      y = MD_mean_regions,
      colour = sex
    ),
    alpha = 0.40,
    size = 1.8
  ) +
  
  # Trayectorias cuadráticas del mismo LMM
  geom_line(
    data = pred_data,
    aes(
      x = pnd,
      y = pred,
      colour = sex,
      group = sex
    ),
    linewidth = 1.3
  ) +
  
  # Mínimos de las trayectorias
  geom_point(
    data = minimos_plot,
    aes(
      x = PND_min,
      y = MD_min,
      colour = sex
    ),
    size = 3.5
  ) +
  
  # Líneas verticales en los mínimos
  geom_vline(
    data = minimos_plot,
    aes(
      xintercept = PND_min,
      colour = sex
    ),
    linetype = "dashed",
    linewidth = 0.7,
    alpha = 0.65
  ) +
  
  # Etiquetas de los mínimos
  geom_text(
    data = minimos_plot,
    aes(
      x = PND_min,
      y = MD_min,
      label = paste0(
        "PND = ", round(PND_min, 1),
        "\nMD = ", round(MD_min, 4)
      ),
      colour = sex
    ),
    vjust = -1,
    size = 3.3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = colores_sex,
    breaks = niveles_sex
  ) +
  
  scale_fill_manual(
    values = colores_sex,
    breaks = niveles_sex
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    title = "Quadratic mixed-model trajectories",
    subtitle = paste0(
      "Mean MD across regions with significant Age × Sex ",
      "interaction after FDR (n = ",
      length(regiones_sig),
      " regions)"
    ),
    x = "Postnatal day (PND)",
    y = "Mean MD across significant regions",
    colour = "Sex",
    fill = "Sex"
  )

# ------------------------------------------------------------
# 7. MOSTRAR EL GRÁFICO
# ------------------------------------------------------------

print(p_media)

# ------------------------------------------------------------
# 8. GUARDAR PDF Y PNG
# ------------------------------------------------------------

carpeta_salida <- file.path(
  "LMM_GM_QUADRATIC_BIC",
  "Figures"
)

dir.create(
  carpeta_salida,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = file.path(
    carpeta_salida,
    "Mean_AgeSex_significant_regions_LMM.pdf"
  ),
  plot = p_media,
  width = 8,
  height = 5.5
)

ggsave(
  filename = file.path(
    carpeta_salida,
    "Mean_AgeSex_significant_regions_LMM.png"
  ),
  plot = p_media,
  width = 8,
  height = 5.5,
  dpi = 300
)

# Confirmación
cat(
  "\nGráfico guardado en:\n",
  normalizePath(
    carpeta_salida,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)

