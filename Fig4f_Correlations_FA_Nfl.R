# ============================================================
# LIBRERÍAS
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)


# ============================================================
# 1. CARGAR DATOS DE FA
# ============================================================

fa <- read_excel(
  "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/FA/NoCleaned_FA_without_outliers/wm_data_FA_without_outliers.xlsx"
)


# ============================================================
# 2. CARGAR DATOS DE NFL
# ============================================================

file <- "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/Alamar_results/data/NULISAseq_NPQ_PP2025_0236_GLIA-MRI.xlsx"

nfl <- read_excel(
  file,
  sheet = "RQ NPQ Values"
)

metadata <- read_excel(
  "D:/Ageing/BLOOD BIOMARKERS/NULISA/Nulisa_analysis/muestras_ensayo.xlsx",
  sheet = 2
)


# ============================================================
# 3. PREPARAR DATOS DE NFL
# ============================================================

df_final <- nfl %>%
  dplyr::select(
    SampleName,
    Target,
    NPQ,
    LOD
  ) %>%
  left_join(
    metadata,
    by = "SampleName"
  ) %>%
  filter(
    !startsWith(AnimalID, "O"),
    !SampleName %in% c(
      "B280222",
      "A230522",
      "B131222",
      "A300123"
    )
  ) %>%
  mutate(
    NPQ = as.numeric(NPQ),
    LOD = as.numeric(LOD),
    detected = NPQ > LOD,
    Sex = as.factor(Sex),
    Timepoint = as.numeric(Timepoint),
    AnimalID = as.factor(AnimalID),
    Protein = Target
  )


# ============================================================
# 4. FILTRO DE PROTEÍNAS SEGÚN LOD
# ============================================================

min_group <- df_final %>%
  group_by(Sex) %>%
  summarise(
    n = n_distinct(SampleName),
    .groups = "drop"
  ) %>%
  summarise(
    min(n)
  ) %>%
  pull()


proteins_keep <- df_final %>%
  group_by(Protein) %>%
  summarise(
    det = sum(detected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    det >= min_group
  ) %>%
  pull(Protein)


df_final <- df_final %>%
  filter(
    Protein %in% proteins_keep
  )


# ============================================================
# 5. PREPARAR EDADES
# ============================================================
#
# FA:
#   Se utilizarán TODAS las edades presentes en fa$pnd.
#
# NfL:
#   Se utilizarán únicamente las edades disponibles:
#
#   161, 280, 364, 546, 609
#
# ============================================================


# ------------------------------------------------------------
# Convertir PND de FA a numérico
# ------------------------------------------------------------

fa <- fa %>%
  mutate(
    pnd = as.numeric(pnd)
  )


# ------------------------------------------------------------
# TODAS LAS EDADES PRESENTES EN FA
# ------------------------------------------------------------

edades_FA <- sort(
  unique(
    fa$pnd[
      !is.na(fa$pnd)
    ]
  )
)


# ------------------------------------------------------------
# EDADES DISPONIBLES EN NFL
# ------------------------------------------------------------

edades_NFL <- c(
  161,
  280,
  364,
  546,
  609
)


# ------------------------------------------------------------
# COMPROBAR EDADES
# ------------------------------------------------------------

print(
  edades_FA
)

print(
  edades_NFL
)


# ============================================================
# 6. PREPARAR NFL
# ============================================================

NFL_todos <- df_final %>%
  filter(
    Protein == "Nefl",
    Timepoint %in% edades_NFL
  ) %>%
  dplyr::select(
    AnimalID,
    SampleName,
    Sex,
    Timepoint,
    NPQ
  )


# ============================================================
# 7. COMPROBAR QUE SOLO HAY UNA MUESTRA POR RATA Y EDAD
# ============================================================

NFL_todos %>%
  count(
    AnimalID,
    Timepoint
  ) %>%
  count(
    n,
    name = "number_of_combinations"
  )


# ============================================================
# 8. COMPROBAR CUÁNTOS ANIMALES HAY EN CADA EDAD DE NFL
# ============================================================

NFL_todos %>%
  count(
    Timepoint
  )


# ============================================================
# 9. CREAR MATRIZ DE CORRELACIONES
#
# FILAS:
#   TODAS las edades de FA presentes en fa$pnd
#
# COLUMNAS:
#   Edades disponibles de NfL
#
# FA:
#   Media de todas las regiones por rata y edad
#
# NFL:
#   Un único valor por rata y edad
#
# Correlación:
#   Spearman
#
# ============================================================

resultados_matrix <- data.frame()


for (edad_FA in edades_FA) {
  
  
  # ==========================================================
  # 9A. FA PARA ESA EDAD
  # ==========================================================
  
  FA_edad <- fa %>%
    dplyr::select(
      experiment,
      rat_code,
      pnd,
      region,
      sex,
      FA
    ) %>%
    filter(
      pnd == edad_FA
    ) %>%
    
    # --------------------------------------------------------
  # MEDIA DE TODAS LAS REGIONES POR RATA
  # --------------------------------------------------------
  
  group_by(
    rat_code
  ) %>%
    
    summarise(
      
      FA_global = ifelse(
        all(is.na(FA)),
        NA_real_,
        mean(
          FA,
          na.rm = TRUE
        )
      ),
      
      sex = first(sex),
      
      n_regions = sum(
        !is.na(FA)
      ),
      
      .groups = "drop"
    )
  
  
  # ==========================================================
  # 9B. COMPROBAR CUÁNTAS REGIONES TIENE CADA RATA
  # ==========================================================
  
  print(
    paste(
      "PND FA:",
      edad_FA
    )
  )
  
  print(
    table(
      FA_edad$n_regions
    )
  )
  
  
  # ==========================================================
  # 9C. RECORRER CADA EDAD DE NFL
  # ==========================================================
  
  for (edad_NFL in edades_NFL) {
    
    
    # --------------------------------------------------------
    # NFL DE ESA EDAD
    # --------------------------------------------------------
    
    NFL_edad <- NFL_todos %>%
      filter(
        Timepoint == edad_NFL
      ) %>%
      dplyr::select(
        AnimalID,
        NPQ
      ) %>%
      filter(
        !is.na(NPQ)
      )
    
    
    # --------------------------------------------------------
    # UNIR SOLO ANIMALES COMUNES
    # --------------------------------------------------------
    
    datos <- FA_edad %>%
      inner_join(
        NFL_edad,
        by = c(
          "rat_code" = "AnimalID"
        )
      ) %>%
      filter(
        !is.na(FA_global),
        !is.na(NPQ)
      )
    
    
    # --------------------------------------------------------
    # CORRELACIÓN SPEARMAN
    # --------------------------------------------------------
    
    if (
      nrow(datos) >= 3 &&
      length(unique(datos$FA_global)) >= 3 &&
      length(unique(datos$NPQ)) >= 3
    ) {
      
      
      test <- cor.test(
        datos$FA_global,
        datos$NPQ,
        method = "spearman",
        exact = FALSE
      )
      
      
      # ------------------------------------------------------
      # GUARDAR RESULTADO
      # ------------------------------------------------------
      
      resultados_matrix <- bind_rows(
        
        resultados_matrix,
        
        tibble(
          FA_age = edad_FA,
          NFL_age = edad_NFL,
          rho = unname(
            test$estimate
          ),
          p = test$p.value,
          n = nrow(datos)
        )
      )
      
      
    } else {
      
      
      # ------------------------------------------------------
      # NO HAY SUFICIENTES DATOS
      # ------------------------------------------------------
      
      resultados_matrix <- bind_rows(
        
        resultados_matrix,
        
        tibble(
          FA_age = edad_FA,
          NFL_age = edad_NFL,
          rho = NA_real_,
          p = NA_real_,
          n = nrow(datos)
        )
      )
      
    }
  }
}


# ============================================================
# 10. CORRECCIÓN FDR
# ============================================================

resultados_matrix <- resultados_matrix %>%
  mutate(
    FDR = p.adjust(
      p,
      method = "fdr"
    )
  )


# ============================================================
# 11. VER TODOS LOS RESULTADOS
# ============================================================

print(
  resultados_matrix
)
write_xlsx(as.data.frame(resultados_matrix), "resultados_matrix_MD_Nfl.xlsx")


# ============================================================
# 12. MATRIZ DE RHO
# ============================================================

matrix_rho <- resultados_matrix %>%
  dplyr::select(
    FA_age,
    NFL_age,
    rho
  ) %>%
  pivot_wider(
    names_from = NFL_age,
    values_from = rho
  ) %>%
  arrange(
    FA_age
  )


print(
  matrix_rho
)


# ============================================================
# 13. MATRIZ DE P
# ============================================================

matrix_p <- resultados_matrix %>%
  dplyr::select(
    FA_age,
    NFL_age,
    p
  ) %>%
  pivot_wider(
    names_from = NFL_age,
    values_from = p
  ) %>%
  arrange(
    FA_age
  )


print(
  matrix_p
)


# ============================================================
# 14. MATRIZ DE FDR
# ============================================================

matrix_FDR <- resultados_matrix %>%
  dplyr::select(
    FA_age,
    NFL_age,
    FDR
  ) %>%
  pivot_wider(
    names_from = NFL_age,
    values_from = FDR
  ) %>%
  arrange(
    FA_age
  )


print(
  matrix_FDR
)


# ============================================================
# 15. SIGNIFICATIVOS SIN CORRECCIÓN
#
# p < 0.05
#
# ============================================================

resultados_significativos <- resultados_matrix %>%
  filter(
    !is.na(p),
    p < 0.05
  ) %>%
  arrange(
    p
  )


print(
  resultados_significativos
)


# ============================================================
# 16. HEATMAP DE SPEARMAN RHO
# ============================================================

p_rho <- ggplot(
  resultados_matrix,
  
  aes(
    x = factor(
      NFL_age,
      levels = edades_NFL
    ),
    
    y = factor(
      FA_age,
      levels = rev(edades_FA)
    ),
    
    fill = rho
  )
) +
  
  geom_tile(
    color = "white"
  ) +
  
  geom_text(
    aes(
      label = ifelse(
        is.na(rho),
        "",
        sprintf(
          "%.2f",
          rho
        )
      )
    ),
    
    size = 3.5
  ) +
  
  scale_x_discrete(
    labels = edades_NFL
  ) +
  
  scale_y_discrete(
    labels = rev(edades_FA)
  ) +
  
  scale_fill_gradient2(
    low = "#08519C",
    mid = "white",
    high = "#CB181D",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  
  labs(
    title = "Mean FA × NfL across ages",
    x = "NfL age (PND)",
    y = "FA age (PND)",
    fill = "Spearman rho"
  ) +
  
  theme_classic(
    base_size = 14
  )


print(
  p_rho
)

# ============================================================
# 17. HEATMAP DE RHO + SIGNIFICACIÓN SIN CORREGIR
#
# * = p < 0.05
#
# IMPORTANTE:
# Aquí NO se utiliza FDR para marcar significación.
#
# ============================================================

p_rho_uncorrected <- ggplot(
  resultados_matrix,
  
  aes(
    x = factor(
      NFL_age,
      levels = edades_NFL
    ),
    
    y = factor(
      FA_age,
      levels = rev(edades_FA)
    ),
    
    fill = rho
  )
) +
  
  geom_tile(
    color = "white"
  ) +
  
  geom_text(
    aes(
      label = case_when(
        
        is.na(rho) ~ "",
        
        !is.na(p) & p < 0.05 ~ paste0(
          sprintf(
            "%.2f",
            rho
          ),
          "*"
        ),
        
        TRUE ~ sprintf(
          "%.2f",
          rho
        )
      )
    ),
    
    size = 3.5
  ) +
  
  scale_x_discrete(
    labels = edades_NFL
  ) +
  
  scale_y_discrete(
    labels = rev(edades_FA)
  ) +
  
  scale_fill_gradient2(
    low = "#08519C",
    mid = "white",
    high = "#CB181D",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  
  labs(
    title = "Mean FA × NfL across ages",
    subtitle = "* p < 0.05, uncorrected",
    x = "NfL age (PND)",
    y = "FA age (PND)",
    fill = "Spearman rho"
  ) +
  
  theme_classic(
    base_size = 14
  )


print(
  p_rho_uncorrected
)

# ============================================================
# HEATMAP DE RHO + FDR
#
# * = FDR < 0.05
# ============================================================

p_rho_FDR <- ggplot(
  resultados_matrix,
  
  aes(
    x = factor(
      NFL_age,
      levels = edades_NFL
    ),
    
    y = factor(
      FA_age,
      levels = rev(edades_FA)
    ),
    
    fill = rho
  )
) +
  
  geom_tile(
    color = "white"
  ) +
  
  geom_text(
    aes(
      label = case_when(
        
        is.na(rho) ~ "",
        
        !is.na(FDR) & FDR < 0.05 ~ paste0(
          sprintf("%.2f", rho),
          "*"
        ),
        
        TRUE ~ sprintf(
          "%.2f",
          rho
        )
      )
    ),
    size = 3.5
  ) +
  
  scale_x_discrete(
    labels = edades_NFL
  ) +
  
  scale_y_discrete(
    labels = rev(edades_FA)
  ) +
  
  scale_fill_gradient2(
    low = "#08519C",
    mid = "white",
    high = "#CB181D",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  
  labs(
    title = "Mean FA × NfL across ages",
    subtitle = "* FDR < 0.05",
    x = "NfL age (PND)",
    y = "FA age (PND)",
    fill = "Spearman rho"
  ) +
  
  theme_classic(
    base_size = 14
  )


print(p_rho_FDR)

