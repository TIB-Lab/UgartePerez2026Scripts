# ============================================================
# MD - PND210 vs PND546
# FEMALE Y MALE POR SEPARADO
# ============================================================

rm(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. CARGAR DATOS
# ============================================================

ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/MD/CleanedMD_without_Outliers/gm_data_MD_without_outliers.xlsx"

md_data <- read_excel(ruta_archivo)

md_data$pnd <- as.numeric(md_data$pnd)

md_data <- md_data %>%
  mutate(
    pnd = if_else(is.na(pnd), 731, pnd)
  )

md_data$region <- as.factor(md_data$region)
md_data$rat_code <- as.factor(md_data$rat_code)
md_data$sex <- as.factor(md_data$sex)


# ============================================================
# 2. REGIONES
# ============================================================

regiones_anteriores <- c(
  "56_GM_m_crtx",
  "35_GM_orbital_cr",
  "68_GM_prelim_crt",
  "67_GM_infralim_c"
)

regiones_intermedias <- c(
  "44_GM_p_a_crtx",
  "37_GM_ss_crtx"
)

regiones_grafico <- c(
  regiones_anteriores,
  regiones_intermedias
)


# Comprobar
setdiff(
  regiones_grafico,
  unique(as.character(md_data$region))
)

# ============================================================
# 3. SELECCIONAR EDADES
# ============================================================

datos_plot <- md_data %>%
  filter(
    region %in% regiones_grafico,
    pnd %in% c(210, 546),
    !is.na(MD),
    !is.na(sex)
  )
# ============================================================
# 4. DEFINIR ÁREA CORTICAL
# ============================================================

datos_plot <- datos_plot %>%
  mutate(
    cortical_area = case_when(
      
      region %in% regiones_anteriores ~
        "anterior cortical areas",
      
      region %in% regiones_intermedias ~
        "intermediate cortical areas"
    )
  )


# ============================================================
# 5. MEDIA DE LAS REGIONES POR RATA
# ============================================================

datos_rat <- datos_plot %>%
  group_by(
    rat_code,
    sex,
    pnd,
    cortical_area
  ) %>%
  summarise(
    MD = mean(MD, na.rm = TRUE),
    n_regiones = sum(!is.na(MD)),
    .groups = "drop"
  )
# ============================================================
# 6. VARIABLES PARA EL GRÁFICO
# ============================================================

datos_rat <- datos_rat %>%
  mutate(
    
    sex = factor(
      sex,
      levels = c("Female", "Male")
    ),
    
    pnd = factor(
      pnd,
      levels = c(210, 546)
    ),
    
    cortical_area = factor(
      cortical_area,
      levels = c(
        "anterior cortical areas",
        "intermediate cortical areas"
      )
    )
  )

# ============================================================
# 9. MEDIA + SEM
# ============================================================

summary_data <- datos_rat %>%
  group_by(
    cortical_area,
    sex,
    pnd
  ) %>%
  summarise(
    
    mean_MD = mean(MD, na.rm = TRUE),
    
    sd_MD = sd(MD, na.rm = TRUE),
    
    n = sum(!is.na(MD)),
    
    se_MD = sd_MD / sqrt(n),
    
    .groups = "drop"
  )
# ============================================================
# 10. POSICIONES X
# ============================================================

summary_data <- summary_data %>%
  mutate(
    
    x = case_when(
      
      cortical_area == "anterior cortical areas" &
        sex == "Female" &
        pnd == 210 ~ 0.85,
      
      cortical_area == "anterior cortical areas" &
        sex == "Female" &
        pnd == 546 ~ 1.15,
      
      cortical_area == "anterior cortical areas" &
        sex == "Male" &
        pnd == 210 ~ 1.85,
      
      cortical_area == "anterior cortical areas" &
        sex == "Male" &
        pnd == 546 ~ 2.15,
      
      cortical_area == "intermediate cortical areas" &
        sex == "Female" &
        pnd == 210 ~ 3.85,
      
      cortical_area == "intermediate cortical areas" &
        sex == "Female" &
        pnd == 546 ~ 4.15,
      
      cortical_area == "intermediate cortical areas" &
        sex == "Male" &
        pnd == 210 ~ 4.85,
      
      cortical_area == "intermediate cortical areas" &
        sex == "Male" &
        pnd == 546 ~ 5.15
    )
  )


# Hacer exactamente las mismas posiciones para los puntos
datos_rat <- datos_rat %>%
  mutate(
    
    x = case_when(
      
      cortical_area == "anterior cortical areas" &
        sex == "Female" &
        pnd == 210 ~ 0.85,
      
      cortical_area == "anterior cortical areas" &
        sex == "Female" &
        pnd == 546 ~ 1.15,
      
      cortical_area == "anterior cortical areas" &
        sex == "Male" &
        pnd == 210 ~ 1.85,
      
      cortical_area == "anterior cortical areas" &
        sex == "Male" &
        pnd == 546 ~ 2.15,
      
      cortical_area == "intermediate cortical areas" &
        sex == "Female" &
        pnd == 210 ~ 3.85,
      
      cortical_area == "intermediate cortical areas" &
        sex == "Female" &
        pnd == 546 ~ 4.15,
      
      cortical_area == "intermediate cortical areas" &
        sex == "Male" &
        pnd == 210 ~ 4.85,
      
      cortical_area == "intermediate cortical areas" &
        sex == "Male" &
        pnd == 546 ~ 5.15
    )
  )
# ============================================================
# 11. COLORES
# ============================================================

colores_MD <- c(
  "Female.210" = "#C8B5F5",
  "Female.546" = "#7650C8",
  "Male.210"   = "#8FE28F",
  "Male.546"   = "#35A84A"
)


summary_data <- summary_data %>%
  mutate(
    grupo = paste(
      sex,
      pnd,
      sep = "."
    )
  )

datos_rat <- datos_rat %>%
  mutate(
    grupo = paste(
      sex,
      pnd,
      sep = "."
    )
  )
# ============================================================
# 12. GRÁFICO FINAL
# ============================================================

p_MD <- ggplot() +
  
  # ----------------------------------------------------------
# BARRAS
# ----------------------------------------------------------

geom_col(
  data = summary_data,
  aes(
    x = x,
    y = mean_MD,
    fill = grupo
  ),
  width = 0.25
) +
  
  
  # ----------------------------------------------------------
# ERROR BARS
# ----------------------------------------------------------

geom_errorbar(
  data = summary_data,
  aes(
    x = x,
    ymin = mean_MD - se_MD,
    ymax = mean_MD + se_MD
  ),
  width = 0.08,
  linewidth = 0.6
) +
  
  
  # ----------------------------------------------------------
# PUNTOS INDIVIDUALES
# ----------------------------------------------------------

geom_jitter(
  data = datos_rat,
  aes(
    x = x,
    y = MD,
    colour = grupo
  ),
  width = 0.07,
  height = 0,
  size = 2.2,
  alpha = 0.65
) +
  
  
  # ----------------------------------------------------------
# COLORES
# ----------------------------------------------------------

scale_fill_manual(
  values = colores_MD,
  name = NULL,
  labels = c(
    "Female.210" = "PND210",
    "Female.546" = "PND546",
    "Male.210" = "PND210",
    "Male.546" = "PND546"
  )
) +
  
  scale_colour_manual(
    values = colores_MD,
    guide = "none"
  ) +
  
  
  # ----------------------------------------------------------
# EJE X
# ----------------------------------------------------------

scale_x_continuous(
  breaks = c(
    1,
    2,
    4,
    5
  ),
  
  labels = c(
    "F",
    "M",
    "F",
    "M"
  ),
  
  limits = c(
    0.55,
    5.45
  )
) +
  
  
  # ----------------------------------------------------------
# EJE Y
# ----------------------------------------------------------

labs(
  x = NULL,
  y = "MD"
) +
  
  
  # ----------------------------------------------------------
# TEMA
# ----------------------------------------------------------

theme_classic() +
  
  theme(
    
    axis.title.y = element_text(
      size = 15,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 14,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 11
    ),
    
    legend.position = "top",
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      10,
      15,
      60,
      10
    )
  )
# ============================================================
# 13. ETIQUETAS DE LOS BLOQUES
# ============================================================

y_min <- -0.12 * max(
  summary_data$mean_MD,
  na.rm = TRUE
)

p_MD <- p_MD +
  
  annotate(
    "text",
    x = 1.5,
    y = y_min,
    label = "anterior\ncortical areas",
    size = 4.8,
    fontface = "bold",
    lineheight = 0.9
  ) +
  
  annotate(
    "text",
    x = 4.5,
    y = y_min,
    label = "intermediate\ncortical areas",
    size = 4.8,
    fontface = "bold",
    lineheight = 0.9
  ) +
  
  coord_cartesian(
    clip = "off"
  )

print(p_MD)


# ###MODELO MIXTO ES MÁS CORRECTO PORQUE HAY MISSING VALUES

library(lme4)
library(lmerTest)

modelo <- lmer(
  MD ~ sex * pnd * cortical_area + (1 | rat_code),
  data = datos_rat
)
summary(modelo)
anova(modelo)

library(emmeans)

emm <- emmeans(
  modelo,
  ~ pnd | sex * cortical_area
)

posthoc <- pairs(
  emm,
  adjust = "holm"
)

# Report differences, 95% CI and adjusted p-values 
posthoc_results <- summary( posthoc, infer = c(TRUE, TRUE)) 
posthoc_results


