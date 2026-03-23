#' Non-probability Sample
#'
#' This dataset is a cohort (convenience) sample derived from the third cycle of the
#' National Health and Nutrition Examination Survey (NHANES). It is used in the
#' IPWMethod package to demonstrate and evaluate inverse probability weighting
#' methods (IPWM) by comparing pseudo-weighted estimators with benchmark estimators
#' based on the original NHANES survey design weights.
#'
#' @format A data frame with 674 observations and 15 variables:
#' \describe{
#'   \item{d_agecat}{Age category (factor with 4 levels)}
#'   \item{d_marital}{Marital status (factor with 4 levels)}
#'   \item{d_race}{Race category (factor with 4 levels)}
#'   \item{d_educat}{Education level (factor with 5 levels)}
#'   \item{d_empstat}{Employment status (factor with 2 levels)}
#'   \item{d_bmicat2}{BMI category (factor with 5 levels)}
#'   \item{d_smoking}{Smoking status (factor with 3 levels)}
#'   \item{d_arthritis}{Arthritis status (factor with 2 levels)}
#'   \item{d_bronchitis}{Bronchitis status (factor with 2 levels)}
#'   \item{d_liver_comorb}{Liver comorbidity status (factor with 2 levels)}
#'   \item{d_osteoporosis}{Osteoporosis status (factor with 2 levels)}
#'   \item{d_comorbidity}{General comorbidity indicator (factor with 2 levels)}
#'   \item{psa}{Outcome variable: serum prostate-specific antigen level (numeric)}
#'   \item{weight}{NHANES survey design weight (numeric)}
#' }
#'
#' @usage data(sc)
#'
#' @examples
#' data(sc)
#' head(sc)
"sc"
