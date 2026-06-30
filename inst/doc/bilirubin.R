## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  echo = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 4.5,
  fig.align = "center",
  warning = TRUE,
  message = TRUE
)

# CRAN-Safe CLI formatting
options(cli.num_colors = 1)

## ----load-libraries, message=FALSE, warning=FALSE-----------------------------
library(survival)
library(dplyr)
library(ggplot2)
library(patchwork)
library(knitr)
library(cli)
library(OptSurvCutR)

## ----load-data----------------------------------------------------------------
data(pbc, package = "survival")

# Extract complete cases for our predictor, outcome, and clinical covariates
pbc_clean <- na.omit(pbc[, c("time", "status", "bili", "age", "sex", "edema")])

# Define the all-cause mortality survival event (1 or 2 = endpoint, 0 = censored)
pbc_clean$event <- as.integer(pbc_clean$status %in% c(1, 2))

head(pbc_clean)

## ----find-number-BIC----------------------------------------------------------
num_res <- find_cutpoint_number(
  # ==========================================
  # 1. CORE DATA & COVARIATE INPUTS
  # ==========================================
  data = pbc_clean, # The patient dataset
  predictor = "bili", # The continuous biomarker
  outcome_time = "time", # Survival time column
  outcome_event = "event", # Survival event column
  covariates = c("age", "sex", "edema"), # Confounders to control for

  # ==========================================
  # 2. SEARCH ENGINE SETTINGS
  # ==========================================
  method = "genetic", # The evolutionary search engine
  criterion = "AIC", # Evaluates adjusted model fit
  max_cuts = 5, # Test models ranging from 0 to 5 cuts
  nmin = 0.15, # SAFETY LIMIT: Minimum 15% of patients per group

  # ==========================================
  # 3. ADVANCED GENETIC TUNING
  # ==========================================
  max.generations = NULL, # LIFESPAN: Number of evolutionary cycles to run
  pop.size = NULL, # SEARCH PARTY: Number of random guesses per generation
  boundary.enforcement = 2, # EDGE CONTROL: 2 = Soft boundary configuration
  seed = 123 # REPRODUCIBILITY
)

## ----find-number-BIC-summary--------------------------------------------------
summary(num_res)

## ----loess--------------------------------------------------------------------
plot(num_res) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    colour = "darkgrey",
    linetype = "dashed",
    alpha = 0.5
  )

## ----find-cutpoint------------------------------------------------------------
cut_res <- find_cutpoint(
  # ==========================================
  # 1. CORE DATA & COVARIATE INPUTS
  # ==========================================
  data = pbc_clean,
  predictor = "bili",
  outcome_time = "time",
  outcome_event = "event",
  covariates = c("age", "sex", "edema"), # Adjusted during location discovery

  # ==========================================
  # 2. SEARCH ENGINE SETTINGS
  # ==========================================
  method = "genetic",
  criterion = "logrank",
  num_cuts = 3, # Setting 3 cuts for clinical tier separation
  nmin = 0.15, # SAFETY LIMIT: Minimum patients per risk cohort
  n_perm = 20, # Streamlined permutation cycle

  # ==========================================
  # 3. ADVANCED GENETIC TUNING
  # ==========================================
  max.generations = NULL,
  pop.size = NULL,
  boundary.enforcement = 2,
  seed = 123,
  n_cores = 1
)

## ----find-cutpoint-summary----------------------------------------------------
summary(cut_res)

## ----plot_find-cutpoint-------------------------------------------------------
plot(cut_res, type = "distribution") +
  geom_rug(alpha = 0.5) +
  labs(caption = "KDE interpolation reveals the covariate-adjusted thresholds of bilirubin")

## ----validate-cutpoint--------------------------------------------------------
val_res <- validate_cutpoint(
  # ==========================================
  # 1. VALIDATION INPUTS
  # ==========================================
  cutpoint_result = cut_res,
  num_replicates = 30, # 500+ is the standard for final publication
  n_cores = 1,

  # ==========================================
  # 2. ADVANCED SETTINGS (Must match Step 2)
  # ==========================================
  max.generations = NULL,
  pop.size = NULL,
  boundary.enforcement = 2,
  seed = 123
)

## ----validate-cutpoint-summary------------------------------------------------
summary(val_res)

## ----validate-cutpoint-knitr--------------------------------------------------
knitr::kable(val_res$confidence_intervals, digits = 3, caption = "Covariate-Adjusted Stability Intervals for Bilirubin Risk Thresholds")

## ----validate-cutpoint-plot---------------------------------------------------
plot(val_res) +
  labs(caption = "Density interpolation of bootstrap-derived cut-points under adjustment")

## ----group-composition--------------------------------------------------------
final_dataset <- plot(cut_res, return_data = TRUE)

composition_table <- final_dataset %>%
  group_by(group) %>%
  summarise(
    Bilirubin_Range = paste0(
      round(min(factor), 3), " – ", round(max(factor), 3)
    ),
    N = n(),
    Deaths = sum(event),
    Mean_Age = round(mean(age), 1),
    Edema_Rate = round(mean(edema), 2)
  )

knitr::kable(composition_table, caption = "Composition and Covariate Profiles of Discovered Bilirubin Risk Groups")

## ----plot-forest, fig.width=7, fig.height=3.5---------------------------------
plot(cut_res,
  type = "forest",
  main = "Adjusted Mortality Hazard Ratios Relative to Group 1"
)

## ----plot-diagnostic----------------------------------------------------------
plot(cut_res, type = "diagnostic")

## ----plot_validation----------------------------------------------------------
plot_validation(val_res,
  focus_cuts = c(1, 2),
  main = "Upper Tier Threshold Stability Map"
)

## ----plot-km------------------------------------------------------------------
plot(cut_res,
  type = "outcome",
  title = "PBC Survival by Bilirubin Tier (Covariate Adjusted)",
  xlab = "Follow-up Time (Days)",
  ylab = "Overall Survival Probability",
  legend.title = "Bilirubin Tier"
)

## ----export-data--------------------------------------------------------------
stratified_patients <- plot(cut_res, return_data = TRUE)
head(stratified_patients[, c("time", "event", "factor", "age", "sex", "edema", "group")])

## ----session_info-------------------------------------------------------------
sessionInfo()

