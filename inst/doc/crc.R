## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "+",
  fig.width = 7, fig.height = 4, fig.align = "center",
  warning = FALSE, message = TRUE
)
options(cli.num_colors = 1)
# Suppress cli progress bars (validate_cutpoint()'s "Bootstrapping" bar):
# same rationale as bilirubin.Rmd - see that vignette's setup chunk.
options(cli.progress_show_after = Inf)

## ----load-libraries-----------------------------------------------------------
library(dplyr)
library(survival)
library(ggplot2)
library(knitr)
library(OptSurvCutR)

options(cli.num_colors = 1)

## ----configuration------------------------------------------------------------
covariates_to_adjust_for <- c("age", "sex", "differ")

# One minimum group size, used at every step
NMIN <- 0.15

## ----load-data----------------------------------------------------------------
# library(survival), loaded above, already exposes `colon` as a lazy-loaded
# object. An explicit data("colon", package = "survival") call here would
# emit a spurious "data set 'colon' not found" warning in a clean session -
# the dataset lives inside survival's multi-object cancer.rda, and data()'s
# own name-based lookup does not resolve it the same way library() does.
# Confirmed harmless and unnecessary; see the equivalent fix in Figure_3.R.

analysis_data <- colon %>%
  filter(etype == 1) %>%                      # recurrence records only
  select(patient_id = id, time_days = time, status,
         nodes, all_of(covariates_to_adjust_for)) %>%
  mutate(
    time_months  = time_days / 30.4375,
    status_final = ifelse(time_months > 60, 0, status),   # 5-year censoring
    time_final   = pmin(time_months, 60)
  ) %>%
  select(patient_id, time = time_final, status = status_final,
         nodes, all_of(covariates_to_adjust_for)) %>%
  filter(complete.cases(time, status, nodes,
                        across(all_of(covariates_to_adjust_for))))

nrow(analysis_data)

## ----predictor-table----------------------------------------------------------
table(analysis_data$nodes)

## ----find-number--------------------------------------------------------------
number_result <- find_cutpoint_number(
  data = analysis_data, predictor = "nodes",
  outcome_time = "time", outcome_event = "status",
  covariates = covariates_to_adjust_for,
  method = "genetic", criterion = "BIC",
  max_cuts = 4, nmin = NMIN,
  max.generations = NULL, pop.size = NULL,
  boundary.enforcement = 2, seed = 123
)

summary(number_result)

## ----plot-number--------------------------------------------------------------
plot(number_result)

## ----find-cutpoint------------------------------------------------------------
cutpoint_result <- find_cutpoint(
  data = analysis_data, predictor = "nodes",
  outcome_time = "time", outcome_event = "status",
  covariates = covariates_to_adjust_for,
  num_cuts = number_result$optimal_num_cuts,
  method = "systematic", criterion = "logrank",
  nmin = NMIN,
  n_perm = 100,          # low for build speed; use >= 1000 when reporting
  seed = 123, n_cores = 1
)

summary(cutpoint_result)

## ----plot-surface, fig.width=7, fig.height=4.75-------------------------------
plot(cutpoint_result, type = "surface")

## ----plot-distribution--------------------------------------------------------
plot(cutpoint_result, type = "distribution")

## ----validate-----------------------------------------------------------------
validation_result <- validate_cutpoint(
  cutpoint_result = cutpoint_result,
  num_replicates = 30,      # reduced for build speed; use >= 500 when reporting
  n_cores = 1, seed = 123   # single core: vignette builds cannot use parallel workers
)

summary(validation_result)

## ----plot-validation----------------------------------------------------------
plot(validation_result)

## ----plot-validation-2d, fig.width=7, fig.height=4----------------------------
plot_validation(validation_result, focus_cuts = c(1, 2),
                main = "Cuts 1 and 2 across resamples")

## ----group-composition--------------------------------------------------------
grouped <- plot(cutpoint_result, return_data = TRUE)

grouped %>%
  group_by(group) %>%
  summarise(
    Nodes    = paste0(min(factor), " - ", max(factor)),
    N        = n(),
    Events   = sum(event),
    Mean_Age = round(mean(age), 1),
    Grade    = round(mean(differ), 2)
  ) %>%
  kable(caption = "Composition of the three nodal burden groups")

## ----plot-forest, fig.width=7, fig.height=3.5---------------------------------
plot(cutpoint_result, type = "forest",
     reference_group = "G1",
     main = "Adjusted hazard ratios for nodal burden groups")

## ----plot-residuals-----------------------------------------------------------
plot(cutpoint_result, type = "diagnostic")

## ----plot-landmark, fig.width=7, fig.height=4.5-------------------------------
plot(cutpoint_result, type = "landmark", landmark = 12,
     legend.title = "Nodal burden group")

## ----plot-km------------------------------------------------------------------
plot(cutpoint_result, type = "outcome",
     title = "Five-year recurrence-free survival by nodal burden",
     xlab = "Time (months)", ylab = "Recurrence-free survival",
     legend.title = "Nodal burden group",
     legend.labs = c("Low (G1)", "Intermediate (G2)", "High (G3)"))

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

