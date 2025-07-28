library(installr)
library(metafor)
library(dplyr)
library(tidyr)
library(meta)
library(ggplot2)
library(forestplot)


# Calculate the overall effect size(Random effect model)
model <- tryCatch({
  rma(yi = effet_size, vi = variance, data = data, mods = ~ 1 | study_id,method = "REML")
})


# Define photosynthetic parameters
photosynthetic_indicators <- c("Chla", "Chlb", "Total chl", "Carotenoids","Fv/F0", "Fv/Fm", "NPQ", "qP", "Pn", "Tr1", "Gs", "Ci", "Rubisco activity")

results_list <- list()

for (indicator in photosynthetic_indicators) {
  indicator_data <- data%>%
    filter(grepl(indicator, index, ignore.case = TRUE))
  if (nrow(indicator_data) == 0) next
  if (nrow(indicator_data) >= 1) {
    model <- tryCatch({
      rma(yi = effet_size, vi = variance, data = indicator_data, random = ~ 1 | study_id,method = "REML")
    }, error = function(e) NULL)
    if (!is.null(model)) {
      results_list[[paste(indicator, sep = "_")]] <- data.frame(
        Indicator = indicator,
        k = model$k,
        effect_size = model$beta,
        ci_lb = model$ci.lb,
        ci_ub = model$ci.ub,
        I2 = model$I2
      )
    }
  }
}
warnings()

results_df <- bind_rows(results_list)


# Calculate the effect size between groups
# Define dose groups
data <- data %>%
  mutate(
    dose_group = case_when(
      dose <= 0.002 ~ "<0.002 mg/g",
      dose> 0.002 & dose <= 0.02 ~ "0.002-0.02 mg/g",
      dose > 0.02 ~ ">0.02 mg/g",
      TRUE ~ NA_character_
    ),
    dose_group = factor(
      dose_group, 
      levels = c("<0.002 mg/g", "0.002-0.02 mg/g", ">0.02 mg/g")
    )
  )

# Define photosynthetic parameters
photosynthetic_indicators <- c("Chla", "Chlb", "Total chl", "Carotenoids","Fv/F0", "Fv/Fm", "NPQ", "qP", "Pn", "Tr1", "Gs", "Ci", "Rubisco activity")

results_list <- list()

for (indicator in photosynthetic_indicators) {
  indicator_data <- data %>%
    filter(grepl(indicator, index, ignore.case = TRUE))
  if (nrow(indicator_data) == 0) next
  for (dose_level in levels(indicator_data$dose_group)) {
    group_data <- indicator_data %>% 
      filter(dose_group == dose_level)
    cat("Processing", indicator, "in", dose_level, "group with", nrow(group_data), "rows.\n")
    if (any(is.na(group_data$effect_size)) || any(is.na(group_data$variance))) {
      cat("Missing values found in", indicator, "in", dose_level, "group. Skipping...\n")
      next
    }
    if (nrow(group_data) >= 1) {
      model <- tryCatch({
        rma(yi = effet_size, vi = variance, data = group_data, mods = ~ 1 | study_id,method = "REML")
      }, error = function(e)error = function(e) {
        cat("Error for", indicator, "in", dose_level, "group:", e$message, "\n")
        NULL
      } )
      if (!is.null(model)) {
        results_list[[paste(indicator, dose_level, sep = "_")]] <- data.frame(
          Indicator = indicator,
          Dose_Group = dose_level,
          k = model$k,
          effect_size = model$beta,
          ci_lb = model$ci.lb,
          ci_ub = model$ci.ub,
          I2 = model$I2
        )
      }
    }
  }
}
warnings()

results_df <- bind_rows(results_list)
