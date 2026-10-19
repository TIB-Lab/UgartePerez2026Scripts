# ============================================================
# BRAIN SIZE ~ AGE + SEX + BODY WEIGHT
# Longitudinal mixed-effects model
#
# Brain size:
#   brain_sizes.csv -> brain_size_voxels
#
# Body weight / age / sex:
#   body_weight_Antonio.csv
# ============================================================


# ============================================================
# 1. PAQUETES
# ============================================================

library(tidyverse)
library(lme4)
library(lmerTest)


# ============================================================
# 2. CARGAR LOS DOS ARCHIVOS
# ============================================================

# Archivo con edad, sexo y peso
datos_body <- read.csv(
  "D:/Ageing/Ampliación_paper_reviews/WM_MRI/body_weight_Antonio.csv",
  stringsAsFactors = FALSE
)

# Archivo con brain size
datos_brain <- read.csv(
  "D:/Ageing/Ampliación_paper_reviews/WM_MRI/brain_sizes.csv",
  stringsAsFactors = FALSE
)


# Comprobar estructura
str(datos_body)
str(datos_brain)

head(datos_body)
head(datos_brain)


# ============================================================
# 3. COMPROBAR LAS VARIABLES
# ============================================================

# Variables del archivo body weight
names(datos_body)

# Variables del archivo brain size
names(datos_brain)


# ============================================================
# 4. PREPARAR EL ARCHIVO DE BRAIN SIZE
# ============================================================

# En brain_sizes.csv:
#
# rat_code       = identificador que corresponde a experiment
# brain_size_voxels = tamaño cerebral

datos_brain <- datos_brain %>%
  select(
    rat_code,
    brain_size_voxels
  ) %>%
  mutate(
    rat_code = as.character(rat_code),
    brain_size_voxels = as.numeric(brain_size_voxels)
  )


# Comprobar si hay varios valores de brain size
# para el mismo rat_code

datos_brain %>%
  group_by(rat_code) %>%
  summarise(
    n_brain_size = n_distinct(brain_size_voxels),
    .groups = "drop"
  ) %>%
  filter(n_brain_size > 1)


# ============================================================
# 5. PREPARAR EL ARCHIVO DE BODY WEIGHT
# ============================================================

datos_body <- datos_body %>%
  mutate(
    experiment = as.character(experiment),
    rat_code = as.character(rat_code),
    pnd = as.numeric(pnd),
    sex = as.numeric(sex),
    weight = as.numeric(weight)
  )


# ============================================================
# 6. UNIR LOS DOS ARCHIVOS
# ============================================================

# MUY IMPORTANTE:
#
# experiment (body_weight_Antonio.csv)
#        =
# rat_code (brain_sizes.csv)
#
# Por eso:

datos <- datos_body %>%
  left_join(
    datos_brain,
    by = c("experiment" = "rat_code")
  )


# ============================================================
# 7. COMPROBAR EL RESULTADO DEL JOIN
# ============================================================

head(datos)

names(datos)

summary(datos$brain_size_voxels)


# Cuántas observaciones tienen brain size
sum(is.finite(datos$brain_size_voxels))

# Cuántas NO tienen brain size
sum(!is.finite(datos$brain_size_voxels))


# ============================================================
# 8. COMPROBAR QUE EL JOIN NO HAYA DUPLICADO OBSERVACIONES
# ============================================================

datos %>%
  group_by(experiment) %>%
  summarise(
    n_brain_values = n_distinct(brain_size_voxels),
    .groups = "drop"
  ) %>%
  filter(n_brain_values > 1)


# ============================================================
# 9. ELIMINAR DUPLICADOS TÉCNICOS
# ============================================================

# Como el archivo original contiene múltiples filas idénticas,
# no queremos tratarlas como observaciones independientes.

datos_modelo <- datos %>%
  distinct(
    experiment,
    rat_code,
    pnd,
    sex,
    weight,
    brain_size_voxels
  )
datos_modelo <- datos_modelo %>%
  filter(!is.na(brain_size_voxels))

# Comprobar cuánto hemos reducido
cat(
  "Filas originales:", nrow(datos), "\n",
  "Filas después de eliminar duplicados:", nrow(datos_modelo), "\n"
)


# ============================================================
# 10. PREPARAR VARIABLES PARA EL MODELO
# ============================================================

datos_modelo <- datos_modelo %>%
  mutate(
    
    # Animal
    rat_code = factor(rat_code),
    
    # Sexo
    # 0 = Male
    # 1 = Female
    sex = factor(
      sex,
      levels = c(0, 1),
      labels = c("Male", "Female")
    ),
    
    # Variables numéricas
    pnd = as.numeric(pnd),
    weight = as.numeric(weight),
    brain_size_voxels = as.numeric(brain_size_voxels)
    
  ) %>%
  
  filter(
    is.finite(pnd),
    is.finite(weight),
    is.finite(brain_size_voxels)
  )


# ============================================================
# 11. CENTRAR EDAD Y PESO
# ============================================================

mean_pnd <- mean(
  datos_modelo$pnd,
  na.rm = TRUE
)

mean_weight <- mean(
  datos_modelo$weight,
  na.rm = TRUE
)


datos_modelo <- datos_modelo %>%
  mutate(
    
    # Edad centrada
    pnd_c = pnd - mean_pnd,
    
    # Peso centrado
    weight_c = weight - mean_weight
    
  )


cat(
  "Mean PND =", mean_pnd, "\n",
  "Mean body weight =", mean_weight, "\n"
)


# ============================================================
# 12. COMPROBAR CUÁNTAS EDADES TIENE CADA ANIMAL
# ============================================================

datos_modelo %>%
  group_by(rat_code) %>%
  summarise(
    n_timepoints = n_distinct(pnd),
    ages = paste(
      sort(unique(pnd)),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# ============================================================
# 13. MODELO MIXTO
# ============================================================

modelo <- lmer(
  brain_size_voxels ~
    pnd_c * sex +
    I(pnd_c^2) +
    weight_c +
    (1 | rat_code),
  data = datos_modelo,
  REML = TRUE
)

summary(modelo)
anova(modelo)

library(ggplot2)
library(dplyr)

# Valores de PND para las predicciones
newdata_brain <- expand.grid(
  pnd_c = seq(
    min(datos_modelo$pnd_c, na.rm = TRUE),
    max(datos_modelo$pnd_c, na.rm = TRUE),
    length.out = 200
  ),
  sex = levels(datos_modelo$sex)
)

# Predicciones del modelo
newdata_brain$predicted <- predict(
  modelo_brain,
  newdata = newdata_brain,
  re.form = NA
)

ggplot(datos_modelo,
       aes(x = pnd_c,
           y = brain_size_voxels,
           color = sex)) +
  
  geom_point(alpha = 0.35, size = 1.8) +
  
  geom_line(
    data = newdata_brain,
    aes(y = predicted),
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = c(
      "Male" = "green",
      "Female" = "purple"
    )
  ) +
  
  theme_classic() +
  
  labs(
    x = "Age",
    y = "Brain size (voxels)",
    color = "Sex"
  )
