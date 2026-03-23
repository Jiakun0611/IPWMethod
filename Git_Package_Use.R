## install IPWMethod from GitHub
install.packages("devtools") #(run once only)
library("survey")
devtools::install_github("Jiakun0611/IPWMethod")


############################################################
# IPWMethod package: examples for built-in data and main functions
#
# This script demonstrates:
#   1. load the example datasets
#   2. create survey design objects
#   3. fit one-reference IPW models
#   4. estimate target quantities
#   5. handle missing values
#   6. fit multi-reference models
#   7. use the optional zcol argument
############################################################

## =========================================================
## 0. Load packages
## =========================================================
install.packages("devtools") #(run only at first time)
library(survey)
devtools::install_github("Jiakun0611/IPWMethod")
library(IPWMethod)
## =========================================================
## 1. Load built-in example datasets
## =========================================================
# These datasets are included in the package:
#   - sc                     : convenience sample (complete case) (NHANES 2000-2004)
#   - sc_with_na             : convenience sample with missing values
#   - ref_survey_1_psu       : probability sample with PSU/strata design (NHANES 2004-2010)
#   - ref_survey_1_bootstrap : probability sample with bootstrap replicate weights
#   - ref_survey_2           : second probability sample (NHIS)

data("sc", package = "IPWMethod")
data("sc_with_na", package = "IPWMethod")
data("ref_survey_1_psu", package = "IPWMethod")
data("ref_survey_1_bootstrap", package = "IPWMethod")
data("ref_survey_2", package = "IPWMethod")

# For convenience, assign shorter names
sp1 <- ref_survey_1_psu
sp2 <- ref_survey_2
sp1_bootstrap <- ref_survey_1_bootstrap

## Quick look at the datasets
head(sc)
head(sp1)
head(sp2)
head(sp1_bootstrap[,1:20])

## =========================================================
## 2. Create survey design objects
## =========================================================

# survey design with PSU and strata
des1 <- svydesign(
  ids = ~psu,
  strata = ~strata,
  weights = ~wts_sp1,
  data = sp1,
  nest = TRUE
)

des2 <- svydesign(
  ids = ~psu_sp2,
  strata = ~strata_sp2,
  weights = ~wts_sp2,
  data = sp2,
  nest = TRUE
)

# bootstrap replicate-weight design
des3 <- svrepdesign(
  data = sp1_bootstrap,
  weights = ~wts_sp1,
  repweights = "^bw[0-9]+$",
  type = "bootstrap"
)

## =========================================================
## 3. One-reference build: p_formula not given
## =========================================================
# The 'method' argument is NOT case sensitive.
# Supported one-reference methods include:
#   - "alp"
#   - "clw"
#   - "cali" / "calibration"

## ---- ALP with missing p_formula ----
fit_alp <- IPWM_build(
  data = list(sc, des1),
  method = "alp"
)
summary(fit_alp)

## ---- Cali with missing p_formula ----
fit_cali <- IPWM_build(
  data = list(sc, des2),
  method = "cali"
)
summary(fit_cali)

## ---- CLW with missing p_formula ----
fit_clw <- IPWM_build(
  data = list(sc, des3),
  method = "CLW"
)
summary(fit_clw)

## =========================================================
## 4. Estimation after build
## =========================================================
est1 <- IPWM_estimate(fit_clw, y = "psa")
summary(est1)

## simple output
IPWM_estimate(fit_clw, y = "psa")

## =========================================================
## 5. One-reference build: user provides p_formula
## =========================================================

fit_cali_formula <- IPWM_build(
  data = list(sc, des1),
  method = "cali",
  p_formula = ~ agecat + marital + race + empstat +
    smoking + comorbidity + psa
)
summary(fit_cali_formula)

## =========================================================
## 6. Missing data example
## =========================================================
# the built-in dataset 'sc_with_na' contains missing values.
# argument na.action controls how missing values are handled in the build step.
# argument sc_wname gives the name of the created pseudo-weight column.

# convenience sample with missing values (2 obs with NA in p_formula, 2 obs with NA in y)
head(sc_with_na)

#   na.action = na.exclude
fit_na_exclude <- IPWM_build(
  data = list(sc_with_na, des2),
  method = "ALP",
  p_formula = ~ agecat + marital + race + empstat +
    smoking + comorbidity + height + weight,
  na.action = na.exclude,
  sc_wname = "pseudo_wts"
)
head(fit_na$sc_updated)

# na.action = na.omit
fit_na_omit <- IPWM_build(
  data = list(sc_with_na, des2),
  method = "ALP",
  p_formula = ~ agecat + marital + race + empstat +
    smoking + comorbidity + height + weight,
  na.action = na.omit,
  sc_wname = "pseudo_wts"
)
head(fit_na$sc_updated)


## =========================================================
## 7. Multi-reference build
## =========================================================
# When more than one reference survey is supplied,
# IPWM_build() fits the multi-reference version.

## ---- 7.1 Multi-reference with all default settings ----
fit_multi_1 <- IPWM_build(
  data = list(sc, des1, des2)
)
summary(fit_multi_1)

est_multi_1 <- IPWM_estimate(fit_multi_1, y = "psa")
summary(est_multi_1)

## ---- 7.2 Multi-reference with user-specified sp_order ----
fit_multi_2 <- IPWM_build(
  data = list(sc, des1, des2),
  precali = TRUE,
  sp_order = "given"
)
summary(fit_multi_2)

## ---- 7.3 Multi-reference with named surveys ----
fit_multi_3 <- IPWM_build(
  data = list(sc, LFS = des1, CTADS = des2)
)
summary(fit_multi_3)

## ---- 7.4 Multi-reference with precali = FALSE ----
fit_multi_4 <- IPWM_build(
  data = list(sc, des1, des2),
  p_formula = list(
    ~ agecat + marital,
    ~ empstat + smoking + comorbidity
  ),
  precali = FALSE
)
summary(fit_multi_4)

## ---- 7.5 Multi-reference using a replicate-weight design + NA handling ----
fit_multi_5 <- IPWM_build(
  data = list(sc_with_na, des3, des2),
  p_formula = list(
    ~ psa,
    ~ agecat + marital + race + empstat +
      smoking + comorbidity + height + weight
  ),
  precali = TRUE,
  sp_order = "size",
  na.action = na.exclude,
  sc_wname = "wts"
)
head(fit_multi_5$sc_updated)
summary(fit_multi_5)

## =========================================================
## 8. Example of zcol
## =========================================================
# Here we show that a binary variable can work whether as a factor
# or as a numeric 0/1 variable in the build stage.

sc_zcol <- sc
sc_zcol$empstat <- as.numeric(sc_zcol$empstat) - 1
str(sc_zcol)  # see that empstat is numeric 0/1 in sc_zcol

fit_zcol_num <- IPWM_build(
  data = list(sc_zcol, des2),
  method = "alp"
)

fit_zcol_factor <- IPWM_build(
  data = list(sc, des2),
  method = "alp"
)

all(fit_zcol_factor$pseudo_weights == fit_zcol_num$pseudo_weights)  # TRUE

