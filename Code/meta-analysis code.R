library(metafor)
library(dplyr)
library(tidyr)
library(ggplot2)
library(multcomp)
library(stringr)
library(openxlsx)
library(gt) 


data <- read.csv('meta-analysis dataset.csv', fileEncoding = "GBK")
data

# publication bias
res_simple <- rma(yi = effect_size, vi = variance, data = data, method = "REML")
summary(res_simple)   
funnel(res_simple, main="Funnel plot")
regtest(res_simple) 
# Trim-and-fill
tf <- trimfill(res_simple,main="Trim-and-fill")
summary(tf)
funnel(tf,main = "Trim-and-fill")


egger_subgroup <- function(dat, group_col) {
  
  groups <- unique(dat[[group_col]])
  
  results <- data.frame(
    Response_variables = character(),
    Z = numeric(),
    p = numeric(),
    df = integer(),
    stringsAsFactors = FALSE
  )
  
  for (g in groups) {
    sub_dat <- dat %>% filter(.data[[group_col]] == g)
    if (nrow(sub_dat) < 3) next
    res <- try(
      rma(yi = effect_size,
          vi =  variance, data = sub_dat, method = "REML"),
      silent = TRUE
    )
    if (inherits(res, "try-error")) next
    egger <- try(
      regtest(res, model = "rma"),
      silent = TRUE
    )
    if (inherits(egger, "try-error")) next
    results <- rbind(
      results,
      data.frame(
        Response_variables = g,
        Z = round(egger$zval, 3),
        p = round(egger$pval, 4),
        df = res$k - 2
      )
    )
  }
  
  return(results)
}

egger_table <- egger_subgroup(
  dat,
  group_col = "Photosynthetic.and.biochemical.indexes"
)


# multilevel randon effect model
model <- tryCatch({
  rma.mv(yi = effect_size,
         V = variance,
         random = list(
           ~ 1 | study_id       
         ),
         method = "REML",
         data = data)
})
summary(model) 


# subgroup analysis
photosynthetic_indicators <- c("Light Reaction", "Calvin Cycle", "Simultaneous Effect")
data$indicator <- factor(data$index group)

model_reg <- rma.mv(
  yi = effect_size,
  V = variance,
  mods = ~ indicator - 1,   
  random = list(
    ~1 | study_id
  ),
  method = "REML",
  data = data
)

summary(model_reg)


data <- data %>%
  mutate(
    dose_group = case_when(
      nanoplastics.dose.mg.g. < 0.002 ~ "<0.002 mg/g",
      nanoplastics.dose.mg.g. >= 0.002 & nanoplastics.dose.mg.g. <=0.02 ~ "0.002-0.02 mg/g",
      nanoplastics.dose.mg.g. > 0.02~ ">0.02 mg/g",
      TRUE ~ NA_character_
    ),
    dose_group = factor(dose_group,
                        levels = c("<0.002 mg/g", "0.002-0.02 mg/g", ">0.02 mg/g"))
  )

data$indicator <- factor(data$index group
photosynthetic_indicators <- c("Light Reaction", "Calvin Cycle", "Simultaneous Effect")
data$mod_group <- interaction(data$index, data$dose_group, sep = "###")
data$mod_group <- droplevels(data$mod_group)

model_reg <- rma.mv(
  yi = effect_size,
  V = variance,
  mods = ~ mod_group - 1,
  random = list(~1 | study_id),
  method = "REML",
  data = data
)
summary(model_reg)


