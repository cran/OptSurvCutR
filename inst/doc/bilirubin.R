## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, echo = TRUE, comment = "+",
  fig.width = 8, fig.height = 4.5, fig.align = "center",
  warning = TRUE, message = TRUE
)
options(cli.num_colors = 1)
# Suppress cli progress bars (validate_cutpoint()'s "Bootstrapping" bar):
# a live progress bar re-renders on every step during a vignette knit,
# producing hundreds of lines of output in the rendered document rather
# than a live terminal display. Interim measure until every
# validate_cutpoint() call below is updated to pass quiet = TRUE directly.
options(cli.progress_show_after = Inf)

## ----load-libraries, message=FALSE, warning=FALSE-----------------------------
library(survival)
library(dplyr)
library(ggplot2)
library(knitr)
library(OptSurvCutR)

## ----load-data----------------------------------------------------------------
data(pbc, package = "survival")

pbc_clean <- na.omit(pbc[, c("time", "status", "bili", "age", "sex", "edema")])

# Endpoint: transplant-free survival
# status: 0 = censored, 1 = transplant, 2 = death
pbc_clean$event <- as.integer(pbc_clean$status %in% c(1, 2))

nrow(pbc_clean)

## ----find-number--------------------------------------------------------------
num_res <- find_cutpoint_number(
  data = pbc_clean, predictor = "bili",
  outcome_time = "time", outcome_event = "event",
  covariates = c("age", "sex", "edema"),
  method = "genetic", criterion = "AIC",
  max_cuts = 5,
  nmin = 0.15,          # each group holds at least 15% of patients
  max.generations = NULL, pop.size = NULL,
  boundary.enforcement = 2, seed = 123
)

summary(num_res)

## ----criterion-curve----------------------------------------------------------
plot(num_res)

## ----find-cutpoint------------------------------------------------------------
cut_res <- find_cutpoint(
  data = pbc_clean, predictor = "bili",
  outcome_time = "time", outcome_event = "event",
  covariates = c("age", "sex", "edema"),
  method = "genetic", criterion = "logrank",
  num_cuts = num_res$optimal_num_cuts,   # carried from step 1
  nmin = 0.15,
  n_perm = 20,          # low for build speed; use >= 1000 when reporting
  max.generations = NULL, pop.size = NULL,
  boundary.enforcement = 2, seed = 123, n_cores = 1
)

summary(cut_res)

## ----plot-distribution--------------------------------------------------------
plot(cut_res, type = "distribution") +
  geom_rug(alpha = 0.5) +
  labs(caption = "Bilirubin thresholds on the marker distribution")

## ----validate-cutpoint--------------------------------------------------------
val_res <- validate_cutpoint(
  cutpoint_result = cut_res,
  num_replicates = 30,      # reduced for build speed; use >= 500 when reporting
  n_cores = 1,      # single core: vignette builds cannot use parallel workers
  max.generations = NULL, pop.size = NULL,
  boundary.enforcement = 2, seed = 123
)

summary(val_res)

## ----group-composition--------------------------------------------------------
final_dataset <- plot(cut_res, return_data = TRUE)

final_dataset %>%
  group_by(group) %>%
  summarise(
    Bilirubin = paste0(round(min(factor), 2), " – ", round(max(factor), 2)),
    N = n(),
    Events = sum(event),
    Mean_Age = round(mean(age), 1),
    Edema = round(mean(edema), 2)
  ) %>%
  kable(caption = "Composition of the four-group model")

## ----plot-forest, fig.width=7, fig.height=3.5---------------------------------
plot(cut_res, type = "forest",
     main = "Adjusted hazard ratios relative to group 1")

## ----plot-diagnostic----------------------------------------------------------
plot(cut_res, type = "diagnostic")

## ----plot-km------------------------------------------------------------------
plot(cut_res, type = "outcome",
     title = "Transplant-free survival by bilirubin group",
     xlab = "Follow-up (days)", ylab = "Transplant-free survival",
     legend.title = "Bilirubin group")

## ----plot-validation-2d-12----------------------------------------------------
plot_validation(val_res, focus_cuts = c(1, 2),
                main = "Cuts 1 and 2 across resamples")

## ----plot-validation-2d-23----------------------------------------------------
plot_validation(val_res, focus_cuts = c(2, 3),
                main = "Cuts 2 and 3 across resamples")

## ----session_info-------------------------------------------------------------
sessionInfo()

