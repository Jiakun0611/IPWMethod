#' NHANES Reference Survey Sample (Cycles 1, 2, 4, and 5)
#'
#' This dataset is a probability-based reference sample constructed from four
#' cycles (1, 2, 4, and 5) of the National Health and Nutrition Examination Survey
#' (NHANES). It serves as one of the reference surveys for estimating
#' pseudo-weights in the cohort sample \code{sc} within the IPWMethod package.
#' The known NHANES survey design weights in this dataset provide population-level
#' benchmarks against which the cohort sample is calibrated to reduce and correct
#' selection bias.
#'
#' @format A data frame with 2872 observations and 14 variables:
#' \describe{
#'   \item{d_agecat}{Age category (factor with 4 levels)}
#'   \item{d_marital}{Marital status (factor with 4 levels)}
#'   \item{d_race}{Race category (factor with 5 levels)}
#'   \item{d_educat}{Education level (factor with 5 levels)}
#'   \item{d_empstat}{Employment status (factor with 2 levels)}
#'   \item{d_bmicat2}{BMI category (factor with 5 levels)}
#'   \item{d_smoking}{Smoking status (factor with 3 levels)}
#'   \item{d_arthritis}{Arthritis status (factor with 2 levels)}
#'   \item{d_bronchitis}{Bronchitis status (factor with 2 levels)}
#'   \item{d_liver_comorb}{Liver comorbidity indicator (factor with 2 levels)}
#'   \item{d_osteoporosis}{Osteoporosis status (factor with 2 levels)}
#'   \item{d_comorbidity}{General comorbidity indicator (factor with 2 levels)}
#'   \item{wts}{NHANES survey design weight (numeric)}
#' }
#'
#' @usage data(ref_survey_1)
#'
#' @examples
#' data(ref_survey_1)
#' head(ref_survey_1)
"ref_survey_1"
