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
library(ggeffects)
library(segmented)

# ============================================================
# CARGAR DATOS
# ============================================================

# FA o FR
ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/FA/NoCleaned_FA_without_outliers/wm_data_FA_without_outliers.xlsx"

wm_data <- read_excel(ruta_archivo)

wm_data$pnd <- as.numeric(wm_data$pnd)

wm_data <- wm_data %>%
  mutate(
    pnd = if_else(is.na(pnd), 731, pnd)
  )

# ============================================================
# ELIMINAR REGIÓN QUE SIGUE MEJOR EL MODELO LINEAL
# ============================================================

wm_data <- wm_data %>%
  filter(region != "59_WM_fas_retrof")

# Comprobar regiones que quedan
unique(wm_data$region)

# ============================================================
# CALCULAR FA MEDIA DE LAS REGIONES WM
# ============================================================

wm_mean <- wm_data %>%
  group_by(rat_code, sex, pnd) %>%
  summarise(
    FA = mean(FA, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# VISUALIZACIÓN EXPLORATORIA
# ============================================================

ggplot(
  wm_mean,
  aes(
    x = pnd,
    y = FA,
    colour = sex,
    fill = sex
  )
) +
  
  geom_point(
    alpha = 0.6,
    size = 2
  ) +
  
  geom_smooth(
    method = "loess",
    se = TRUE,
    linewidth = 1.5
  ) +
  
  scale_color_manual(
    values = c(
      "Male" = "#77DD77",
      "Female" = "#9370DB"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Male" = "#77DD77",
      "Female" = "#9370DB"
    )
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Postnatal day",
    y = "Fractional anisotropy",
    colour = "Sex",
    fill = "Sex"
  )

# ============================================================
# MODELO SEGMENTED
# ============================================================

# Separar por sexo
female <- wm_mean %>%
  filter(sex == "Female")

male <- wm_mean %>%
  filter(sex == "Male")

# ============================================================
# MODELOS LINEALES
# ============================================================

mod_f <- lm(
  FA ~ pnd,
  data = female
)

mod_m <- lm(
  FA ~ pnd,
  data = male
)

# ============================================================
# MODELOS SEGMENTADOS
# ============================================================

seg_f <- segmented(
  mod_f,
  seg.Z = ~ pnd,
  psi = 200
)

seg_m <- segmented(
  mod_m,
  seg.Z = ~ pnd,
  psi = 200
)

# ============================================================
# BREAKPOINTS
# ============================================================

bp_female <- seg_f$psi[, "Est."]
bp_male <- seg_m$psi[, "Est."]

cat("Breakpoint Female:", bp_female, "\n")
cat("Breakpoint Male:", bp_male, "\n")

# ============================================================
# MALLA DE PREDICCIÓN
# ============================================================

new_f <- data.frame(
  pnd = seq(
    min(female$pnd),
    max(female$pnd),
    length.out = 200
  )
)

new_m <- data.frame(
  pnd = seq(
    min(male$pnd),
    max(male$pnd),
    length.out = 200
  )
)

# ============================================================
# PREDICCIÓN + IC95%
# ============================================================

pred_f <- predict(
  seg_f,
  newdata = new_f,
  se.fit = TRUE
)

pred_m <- predict(
  seg_m,
  newdata = new_m,
  se.fit = TRUE
)

# ============================================================
# DATAFRAME PREDICCIÓN FEMALE
# ============================================================

pred_f_df <- data.frame(
  pnd = new_f$pnd,
  FA = pred_f$fit,
  lower = pred_f$fit - 1.96 * pred_f$se.fit,
  upper = pred_f$fit + 1.96 * pred_f$se.fit,
  sex = "Female"
)

# ============================================================
# DATAFRAME PREDICCIÓN MALE
# ============================================================

pred_m_df <- data.frame(
  pnd = new_m$pnd,
  FA = pred_m$fit,
  lower = pred_m$fit - 1.96 * pred_m$se.fit,
  upper = pred_m$fit + 1.96 * pred_m$se.fit,
  sex = "Male"
)

# ============================================================
# UNIR PREDICCIONES
# ============================================================

pred <- bind_rows(
  pred_f_df,
  pred_m_df
)

# ============================================================
# FIGURA
# ============================================================

ggplot(
  wm_mean,
  aes(
    x = pnd,
    y = FA,
    colour = sex,
    fill = sex
  )
) +
  
  # Puntos individuales
  geom_point(
    alpha = 0.5,
    size = 2
  ) +
  
  # IC95%
  geom_ribbon(
    data = pred,
    aes(
      x = pnd,
      ymin = lower,
      ymax = upper,
      fill = sex
    ),
    inherit.aes = FALSE,
    alpha = 0.20,
    colour = NA
  ) +
  
  # Modelo segmentado
  geom_line(
    data = pred,
    aes(
      x = pnd,
      y = FA,
      colour = sex
    ),
    inherit.aes = FALSE,
    linewidth = 1.8
  ) +
  
  scale_color_manual(
    values = c(
      "Male" = "#77DD77",
      "Female" = "#9370DB"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Male" = "#77DD77",
      "Female" = "#9370DB"
    )
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Postnatal day",
    y = "Fractional anisotropy",
    colour = "Sex",
    fill = "Sex"
  )


