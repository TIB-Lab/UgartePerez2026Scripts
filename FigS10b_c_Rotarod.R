####MEAN day tests, across age, with weight

rm(list=ls())

# ============================================================
# 1. LEER LOS DOS EXCEL
# ============================================================

df_15 <- read_excel(
  "D:/Ageing/Conducta_Egoa1/BATCH2/rotarod/1.5 años/rats.xlsx"
)

df_20 <- read_excel(
  "D:/Ageing/Conducta_Egoa1/BATCH1&2/2YEARS/rotarod/rats.xlsx"
)


# ============================================================
# 2. PASAR CADA DATASET A FORMATO LARGO
# ============================================================

df_long_15 <- df_15 %>%
  pivot_longer(
    cols = -1,
    names_to = "Column",
    values_to = "Value"
  ) %>%
  rename(Day = ...1) %>%
  mutate(Edad = "1.5")


df_long_20 <- df_20 %>%
  pivot_longer(
    cols = -1,
    names_to = "Column",
    values_to = "Value"
  ) %>%
  rename(Day = ...1) %>%
  mutate(Edad = "2")


# ============================================================
# 3. JUNTAR LOS DOS EXPERIMENTOS
# ============================================================

df_long <- bind_rows(
  df_long_15,
  df_long_20
)


# ============================================================
# 4. SACAR EL PESO DE CADA RATA
# ============================================================

weights <- df_long %>%
  filter(Day == "Weight") %>%
  dplyr::select(Column, Edad, Weight = Value)


# ============================================================
# 5. QUEDARNOS CON DAY1-DAY3
# ============================================================

df_long_days <- df_long %>%
  filter(Day != "Weight") %>%
  left_join(
    weights,
    by = c("Column", "Edad")
  ) %>%
  mutate(
    Day = factor(
      Day,
      levels = c("Day1", "Day2", "Day3")
    ),
    Edad = factor(
      Edad,
      levels = c("1.5", "2")
    ),
    Column = factor(Column)
  )


# ============================================================
# 6. SEXO
# ============================================================

df_long_days <- df_long_days %>%
  mutate(
    Sexo = case_when(
      Column %in% c(
        "Rat4", "Rat5", "Rat6", "Rat7", "Rat8"
      ) ~ "Male",
      
      Column %in% c(
        "Rat9", "Rat10", "Rat11", "Rat12", "Rat13"
      ) ~ "Female",
      
      TRUE ~ NA_character_
    ),
    Sexo = factor(Sexo)
  )


# ============================================================
# 7. CREAR RAT_CODE
# ============================================================

df_long_days <- df_long_days %>%
  mutate(
    rat_code = Column
  )

View(df_long_days)
df_long_days_filtrado <- df_long_days %>%
  filter(!is.na(Day))
# ============================================================
# 8. CALCULAR LA MEDIA DE LOS 3 DIAS POR RATA Y EDAD
# ============================================================

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

df_media <- df_long_days_filtrado %>%
  group_by(rat_code, Column, Sexo, Edad) %>%
  summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    Weight = first(Weight),
    n_dias = sum(!is.na(Value)),
    .groups = "drop"
  )

# Comprobar los datos resumidos
View(df_media)


# ============================================================
# 9. MODELO LINEAL MIXTO CON LA MEDIA DE LOS 3 DIAS
# ============================================================

modelo_media <- lmer(
  Mean_Value ~ Sexo * Edad + Weight +
    (1 | rat_code),
  data = df_media
)

summary(modelo_media)
anova(modelo_media)


# ============================================================
# 10. COMPARAR SEXO Y EDAD
# ============================================================

# Comparación de sexos dentro de cada edad
emmeans(modelo_media, ~ Sexo | Edad)

# Comparación de edades dentro de cada sexo
emmeans(modelo_media, pairwise ~ Edad | Sexo)



#####LEARNING RATE, no weight, across age

rm(list=ls())

# ============================================================
# 1. LEER LOS DOS EXCEL
# ============================================================

df_15 <- read_excel(
  "D:/Ageing/Conducta_Egoa1/BATCH2/rotarod/1.5 años/rats.xlsx"
)

df_20 <- read_excel(
  "D:/Ageing/Conducta_Egoa1/BATCH1&2/2YEARS/rotarod/rats.xlsx"
)


# ============================================================
# 2. PASAR CADA DATASET A FORMATO LARGO
# ============================================================

df_long_15 <- df_15 %>%
  pivot_longer(
    cols = -1,
    names_to = "Column",
    values_to = "Value"
  ) %>%
  rename(Day = ...1) %>%
  mutate(Edad = "1.5")


df_long_20 <- df_20 %>%
  pivot_longer(
    cols = -1,
    names_to = "Column",
    values_to = "Value"
  ) %>%
  rename(Day = ...1) %>%
  mutate(Edad = "2")


# ============================================================
# 3. JUNTAR LOS DOS EXPERIMENTOS
# ============================================================

df_long <- bind_rows(
  df_long_15,
  df_long_20
)


# ============================================================
# 4. SACAR EL PESO DE CADA RATA
# ============================================================

weights <- df_long %>%
  filter(Day == "Weight") %>%
  dplyr::select(Column, Edad, Weight = Value)


# ============================================================
# 5. QUEDARNOS CON DAY1-DAY3
# ============================================================

df_long_days <- df_long %>%
  filter(Day != "Weight") %>%
  left_join(
    weights,
    by = c("Column", "Edad")
  ) %>%
  mutate(
    Day = factor(
      Day,
      levels = c("Day1", "Day2", "Day3")
    ),
    Edad = factor(
      Edad,
      levels = c("1.5", "2")
    ),
    Column = factor(Column)
  )


# ============================================================
# 6. SEXO
# ============================================================

df_long_days <- df_long_days %>%
  mutate(
    Sexo = case_when(
      Column %in% c(
        "Rat4", "Rat5", "Rat6", "Rat7", "Rat8"
      ) ~ "Male",
      
      Column %in% c(
        "Rat9", "Rat10", "Rat11", "Rat12", "Rat13"
      ) ~ "Female",
      
      TRUE ~ NA_character_
    ),
    Sexo = factor(Sexo)
  )


# ============================================================
# 7. CREAR RAT_CODE
# ============================================================

df_long_days <- df_long_days %>%
  mutate(
    rat_code = Column
  )

View(df_long_days)
df_long_days_filtrado <- df_long_days %>%
  filter(!is.na(Day))

modelo <- lmer(
  Value ~ Day * Sexo * Edad +
    (1 | rat_code),
  data = df_long_days_filtrado
)

summary(modelo)
anova(modelo)
library(emmeans)
emmeans(modelo, ~ Sexo * Edad)
emmeans(modelo, pairwise ~ Edad | Sexo)

