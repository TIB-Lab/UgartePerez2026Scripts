library(tidyverse)

# 1. Leer archivos
datos <- read.csv(
  "D:/Ageing/Ampliación_paper_reviews/WM_MRI/body_weight_Antonio.csv"
)

brain <- read.csv(
  "D:/Ageing/Ampliación_paper_reviews/WM_MRI/brain_sizes.csv"
)

# 2. Quedarnos con un único weight por experiment
weights <- datos %>%
  dplyr::select(experiment, weight) %>%
  distinct()

# 3. Unir:
#    experiment (body_weight_Antonio.csv)
#    con rat_code (brain_sizes.csv)
datos_cor <- weights %>%
  inner_join(
    brain,
    by = c("experiment" = "rat_code")
  )

# 4. Comprobar resultado
head(datos_cor)
nrow(datos_cor)

# 5. Asegurar que las variables son numéricas
datos_cor <- datos_cor %>%
  mutate(
    weight = as.numeric(weight),
    brain_size_voxels = as.numeric(brain_size_voxels)
  ) %>%
  filter(
    is.finite(weight),
    is.finite(brain_size_voxels)
  )

# 6. Correlación de Pearson
cor.test(
  datos_cor$weight,
  datos_cor$brain_size_voxels,
  method = "pearson"
)

# 7. Gráfico
ggplot(
  datos_cor,
  aes(x = weight, y = brain_size_voxels)
) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  theme_classic() +
  labs(
    x = "Body weight",
    y = "Brain size (voxels)",
    title = "Relationship between body weight and brain size"
  )
