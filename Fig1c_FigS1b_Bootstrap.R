
rm(list = ls())

library(openxlsx)
library(readxl)
library(dplyr)
library(readr)
library(rsample)
library(modelr)
library(readxl)

#FR
ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/FR/WM/wm_data_FR_without_outliers.xlsx"
#FA
#ruta_archivo <- "D:/Ageing/Modelos_de_regresión_R_STUDIO/results_roi/Exploratory_analysis_Seg_regression/FA/NoCleaned_FA_without_outliers/wm_data_FA_without_outliers.xlsx"

wm_data <- read_xlsx(ruta_archivo)


wm_data <- read_excel(ruta_archivo)
wm_data$pnd <- as.numeric(wm_data$pnd)


wm_data <- wm_data %>%
  mutate(pnd = if_else(is.na(pnd), 731, pnd))

wm_data$region <- as.factor(wm_data$region) 
wm_regiones <- levels(wm_data$region)
wm_regiones <- levels(droplevels(wm_data$region))

#Eliminamos esta región solo en FA porque sigue mejor modelo lineal, en FR no
wm_regiones <- setdiff(
  wm_regiones,
  "59_WM_fas_retrof"
)


####ANIMAL-LEVEL BOOTSTRAP
############################################

library(segmented)
library(dplyr)


#############################
## BOOTSTRAP SETTINGS
#############################

B <- 5000

set.seed(12345)

# Archivo con brain size
datos_brain <- read.csv(
  "D:/Ageing/Ampliación_paper_reviews/WM_MRI/brain_sizes.csv",
  stringsAsFactors = FALSE
)


wm_data <- wm_data %>%
  left_join(
    brain_size,
    by = c("experiment" = "rat_code")
  )

#############################
## FUNCIÓN PARA ESTIMAR BP
#############################

estimate_bp <- function(data, psi0 = 200){
  
  mod <- lm(
    FR ~ pnd + brain_size_voxels,
    data = data
  )
  
  seg <- segmented(
    mod,
    seg.Z = ~pnd,
    psi = psi0
  )
  
  bp <- seg$psi[, "Est."]
  
  return(as.numeric(bp))
}



#############################
## FUNCIÓN BOOTSTRAP REGIÓN
#############################

bootstrap_region <- function(datos_region,
                             B=5000){
  
  
  ids_m <- unique(
    datos_region$rat_code[
      datos_region$sex=="Male"
    ])
  
  
  ids_f <- unique(
    datos_region$rat_code[
      datos_region$sex=="Female"
    ])
  
  
  results <- vector("list", B)
  
  
  for(b in 1:B){
    
    
    #-------------------------
    # Resamplear ratas completas
    #-------------------------
    
    boot_ids_m <- sample(
      ids_m,
      size=length(ids_m),
      replace=TRUE
    )
    
    
    boot_ids_f <- sample(
      ids_f,
      size=length(ids_f),
      replace=TRUE
    )
    
    
    boot_m <- bind_rows(
      
      lapply(seq_along(boot_ids_m),
             
             function(k){
               
               d <- datos_region %>%
                 filter(
                   rat_code==boot_ids_m[k],
                   sex=="Male"
                 )
               
               d$boot_rat <- paste0(
                 "M_",k
               )
               
               d
               
             })
    )
    
    
    
    boot_f <- bind_rows(
      
      lapply(seq_along(boot_ids_f),
             
             function(k){
               
               d <- datos_region %>%
                 filter(
                   rat_code==boot_ids_f[k],
                   sex=="Female"
                 )
               
               d$boot_rat <- paste0(
                 "F_",k
               )
               
               d
               
             })
    )
    
    
    
    boot_data <- bind_rows(
      boot_m,
      boot_f
    )
    
    
    
    #-------------------------
    # Ajustes segmentados
    #-------------------------
    
    res <- tryCatch({
      
      
      bp_m <- estimate_bp(
        boot_data %>%
          filter(sex=="Male"),
        psi0=200
      )
      
      
      bp_f <- estimate_bp(
        boot_data %>%
          filter(sex=="Female"),
        psi0=200
      )
      
      
      data.frame(
        
        iter=b,
        
        bp_male=bp_m,
        
        bp_female=bp_f,
        
        delta=bp_m-bp_f
        
      )
      
      
    },
    
    error = function(e){
      message("ERROR en iteración ", b, ": ", e$message)
      NULL
    })
    
    
    results[[b]] <- res
    
    
    if(b %% 500 == 0){
      cat("Bootstrap:", b, "/", B, "\n")
    }
    
    
  }
  
  
  bind_rows(results)
  
}



#############################
## EJECUTAR TODAS LAS REGIONES
#############################

bootstrap_results_list <- list()


for(reg in wm_regiones){
  
  
  cat("\n====================\n")
  cat("Region:", reg, "\n")
  cat("====================\n")
  
  
  datos_region <- wm_data %>%
    filter(region==reg)
  
  
  tmp <- bootstrap_region(
    datos_region,
    B=B
  )
  
  
  tmp$region <- reg
  
  
  bootstrap_results_list[[reg]] <- tmp
  
}



#############################
## UNIR RESULTADOS
#############################

bootstrap_results <- bind_rows(
  bootstrap_results_list
)

############################################################
## BOOTSTRAP STABILITY
############################################################

B <- 5000

bootstrap_stability <- bootstrap_results %>%
  group_by(region) %>%
  summarise(
    
    ########################################################
    ## 1. CONVERGENCIA
    ########################################################
    
    n_success = sum(
      is.finite(bp_male) &
        is.finite(bp_female) &
        is.finite(delta)
    ),
    
    convergence_rate = n_success / B,
    
    
    ########################################################
    ## 2. ESTABILIDAD DEL BREAKPOINT EN MACHOS
    ########################################################
    
    mean_bp_male = mean(bp_male, na.rm = TRUE),
    
    median_bp_male = median(bp_male, na.rm = TRUE),
    
    sd_bp_male = sd(bp_male, na.rm = TRUE),
    
    IQR_bp_male = IQR(bp_male, na.rm = TRUE),
    
    CI_male_low = quantile(
      bp_male,
      0.025,
      na.rm = TRUE
    ),
    
    CI_male_high = quantile(
      bp_male,
      0.975,
      na.rm = TRUE
    ),
    
    CI_male_width =
      CI_male_high - CI_male_low,
    
    
    ########################################################
    ## 3. ESTABILIDAD DEL BREAKPOINT EN HEMBRAS
    ########################################################
    
    mean_bp_female = mean(bp_female, na.rm = TRUE),
    
    median_bp_female = median(bp_female, na.rm = TRUE),
    
    sd_bp_female = sd(bp_female, na.rm = TRUE),
    
    IQR_bp_female = IQR(bp_female, na.rm = TRUE),
    
    CI_female_low = quantile(
      bp_female,
      0.025,
      na.rm = TRUE
    ),
    
    CI_female_high = quantile(
      bp_female,
      0.975,
      na.rm = TRUE
    ),
    
    CI_female_width =
      CI_female_high - CI_female_low,
    
    
    ########################################################
    ## 4. ESTABILIDAD DE LA DIFERENCIA
    ##    Male - Female
    ########################################################
    
    mean_delta = mean(delta, na.rm = TRUE),
    
    median_delta = median(delta, na.rm = TRUE),
    
    sd_delta = sd(delta, na.rm = TRUE),
    
    IQR_delta = IQR(delta, na.rm = TRUE),
    
    CI_delta_low = quantile(
      delta,
      0.025,
      na.rm = TRUE
    ),
    
    CI_delta_high = quantile(
      delta,
      0.975,
      na.rm = TRUE
    ),
    
    CI_delta_width =
      CI_delta_high - CI_delta_low,
    
    
    ########################################################
    ## 5. ESTABILIDAD DE LA DIRECCIÓN DEL EFECTO
    ########################################################
    
    proportion_negative =
      mean(delta < 0, na.rm = TRUE),
    
    proportion_positive =
      mean(delta > 0, na.rm = TRUE),
    
    ########################################################
    ## 6. PROPORCIÓN DE BOOTSTRAPS QUE CRUZAN 0
    ########################################################
    
    proportion_crossing_zero =
      mean(
        delta <= 0 & delta >= 0,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

#############################
## GUARDAR RESULTADOS CRUDOS
#############################

write.csv(
  bootstrap_results,
  "bootstrap_breakpoints_B5000_raw_delete.csv",
  row.names=FALSE
)
bootstrap_summary <- bootstrap_results %>%
  group_by(region) %>%
  summarise(
    
    mean_bp_male = mean(bp_male),
    mean_bp_female = mean(bp_female),
    
    mean_delta = mean(delta),
    
    CI_low = quantile(delta,0.025),
    CI_high = quantile(delta,0.975),
    
    n_success = n(),
    
    .groups="drop"
  )


write.csv(
  bootstrap_summary,
  "bootstrap_breakpoints_B5000_summary_pvalue_
  FR.csv",
  row.names=FALSE
)

bootstrap_summary


write.csv(
  bootstrap_stability,
  "bootstrap_breakpoints_B5000_stability_FR_WM.csv",
  row.names=FALSE
)

##REPRESENTACION GRAFICA 
####################################################
## BREAKPOINT MEDIO DE TODA LA WM (5000 BOOTSTRAPS)
####################################################

bootstrap_global <- bootstrap_results %>%
  group_by(iter) %>%
  summarise(
    bp_male   = mean(bp_male, na.rm = TRUE),
    bp_female = mean(bp_female, na.rm = TRUE),
    delta     = bp_male - bp_female,
    .groups = "drop"
  )

bootstrap_long <- bootstrap_global %>%
  tidyr::pivot_longer(
    cols = c(bp_male, bp_female),
    names_to = "sex",
    values_to = "breakpoint"
  ) %>%
  mutate(
    sex = recode(
      sex,
      bp_male = "Male",
      bp_female = "Female"
    )
  )

####################################################
## FIGURA BOOTSTRAP 
####################################################

library(ggplot2)

p_boot <- ggplot(
  bootstrap_long,
  aes(
    x = sex,
    y = breakpoint,
    fill = sex
  )
) +
  
  geom_violin(
    trim = FALSE,
    alpha = 0.35,
    colour = NA
  ) +
  
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    fill = "white",
    colour = "black"
  ) +
  
  scale_fill_manual(values = c(
    Male = "#77DD77",
    Female = "#9370DB"
  )) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "",
    y = "Mean WM breakpoint (PND)"
  )

p_boot

ggsave(
  "Bootstrap_mean_breakpoint_WM_FR.pdf",
  p_boot,
  width = 4,
  height = 5
)

####################################################
## BOOTSTRAP SOLO REGIONES CON DIFERENCIA NEGATIVA
####################################################

regiones_seleccionadas <- c(
  "19_WM_cc_1",
  "20_WM_cc_2",
  "5_WM_ic_cerebral"
)


bootstrap_selected <- bootstrap_results %>%
  filter(region %in% regiones_seleccionadas) %>%
  group_by(iter) %>%
  summarise(
    bp_male = mean(bp_male, na.rm = TRUE),
    bp_female = mean(bp_female, na.rm = TRUE),
    delta = bp_male - bp_female,
    .groups = "drop"
  )
bootstrap_selected %>%
  summarise(
    mean_bp_male = mean(bp_male),
    mean_bp_female = mean(bp_female),
    mean_delta = mean(delta),
    CI_low = quantile(delta,0.025),
    CI_high = quantile(delta,0.975)
  )

####################################################
## REPRESENTACIÓN SOLO REGIONES SELECCIONADAS: una figura 
####################################################

bootstrap_long_selected <- bootstrap_selected %>%
  tidyr::pivot_longer(
    cols = c(bp_male, bp_female),
    names_to = "sex",
    values_to = "breakpoint"
  ) %>%
  mutate(
    sex = recode(
      sex,
      bp_male = "Male",
      bp_female = "Female"
    )
  )

p_boot_selected <- ggplot(
  bootstrap_long_selected,
  aes(
    x = sex,
    y = breakpoint,
    fill = sex
  )
) +
  
  geom_boxplot(
    width = 0.35,
    outlier.shape = NA,
    alpha = 0.6
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    fill = "white",
    colour = "black"
  ) +
  
  scale_fill_manual(values = c(
    Male = "#77DD77",
    Female = "#9370DB"
  )) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "",
    y = "Mean breakpoint PND\n(selected WM regions)"
  )


p_boot_selected

##################OTRA FIGURA DEL PANEL

bootstrap_all_regions <- bootstrap_results %>%
  group_by(region) %>%
  summarise(
    mean_delta = mean(delta, na.rm=TRUE),
    CI_low = quantile(delta,0.025, na.rm=TRUE),
    CI_high = quantile(delta,0.975, na.rm=TRUE),
    .groups="drop"
  )


p_all <- ggplot(
  bootstrap_all_regions,
  aes(
    x=reorder(region, mean_delta),
    y=mean_delta
  )
)+
  
  geom_hline(
    yintercept=0,
    linetype="dashed"
  )+
  
  geom_errorbar(
    aes(
      ymin=CI_low,
      ymax=CI_high
    ),
    width=0.2,
    linewidth=0.8
  )+
  
  geom_point(
    size=3
  )+
  
  coord_flip()+
  
  theme_classic(base_size=14)+
  
  labs(
    x="WM region",
    y="Breakpoint difference\n(Male - Female, PND)"
  )

p_all

##bonito
forest_data <- bootstrap_results %>%
  group_by(region) %>%
  reframe(
    estimate = mean(delta, na.rm=TRUE),
    CI_low = as.numeric(quantile(delta,0.025,na.rm=TRUE)),
    CI_high = as.numeric(quantile(delta,0.975,na.rm=TRUE)),
    selected = ifelse(
      region %in% c(
        "5_WM_ic_cerebral_peduncle"
      ),
      "Consistent effect",
      "Other regions"
    )
  ) %>%
  arrange(estimate) %>%
  mutate(
    region=factor(region, levels=unique(region))
  )

ggplot(
  forest_data,
  aes(
    y=region,
    x=estimate,
    colour=selected
  )
)+
  
  # zona sin efecto
  annotate(
    "rect",
    xmin=-10,
    xmax=10,
    ymin=-Inf,
    ymax=Inf,
    fill="grey90",
    alpha=0.5
  )+
  
  geom_vline(
    xintercept=0,
    linetype="dashed",
    colour="grey40"
  )+
  
  geom_errorbarh(
    aes(
      xmin=CI_low,
      xmax=CI_high
    ),
    height=0.25,
    linewidth=1
  )+
  
  geom_point(
    size=4
  )+
  
  
  scale_colour_manual(
    values=c(
      "Other regions"="grey65",
      "Consistent effect"="#D55E00"
    )
  )+
  
  theme_classic(base_size=15)+
  
  theme(
    legend.position="top",
    axis.title.y=element_blank(),
    axis.text.y=element_text(size=12),
    panel.grid.major.y=element_blank()
  )+
  
  labs(
    x="Breakpoint difference (Male - Female, PND)",
    colour=""
  )
