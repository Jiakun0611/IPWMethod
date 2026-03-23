library(survey)


options(survey.lonely.psu = "adjust")


df <- data.frame(
  strata = c(1,1, 2,2, 3),   # lonely PSU
  psu    = c(1,2, 1,2, 1),
  w      = c(1.2, 0.8, 1.5, 1.1, 2.0),
  x      = c(10, 20, 15, 25, 30)
)

df <- data.frame(
  strata = c(1,1,1,1,  2,2,2,2),
  psu    = c(1,1,2,2,  1,1,2,2),  # no lonely psu
  w      = c(1.2,1.1,0.9,1.0,  1.5,1.4,1.3,1.2),
  x      = c(10,12,20,22,  15,17,25,27)
)

des <- svydesign(
  ids = ~psu,
  strata = ~strata,
  weights = ~w,
  data = df,
  nest = TRUE
)

fit <- svytotal(~x, design = des)

help("svyCprod", help_type = "html")
help(package = "survey")

print(fit)
str(fit)
attr(fit, "var")
vcov(fit)
