data <- data.frame(
  psuuu = c(1,1,2,2,3,3),
  strataaa = c(1,1,1,1,2,2),
  weighttt = c(2.1,2.1,1.8,1.8,3.0,3.0),
  income = c(50,60,45,55,70,65),
  age = c(40,50,35,45,60,55)
)

library(survey)

design <- svydesign(
  ids = ~psuuu,
  strata = ~strataaa,
  weights = ~weighttt,
  data = data,
  nest = TRUE
)

class(design)
str(design)
design$variables
weights(design)  # recommended
design$cluster
design$strata
model.frame(design)


