rm(list = setdiff(ls(), c("sc", "sp", "build")))

head(sc); head(sp1); head(sp2)


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

des3 <- svrepdesign(
  data = sp1_bootstrap,
  weights = ~wts_sp1,
  repweights = "^bw[0-9]+$",
  type = "bootstrap"
)

head(des3$variables[,1:20])



# calibration/cali/alp/clw/multi(not case sensitive)
# no formula
fit = IPWM_build(data = list(sc,des1),method ="alp")
summary(fit)

fit = IPWM_build(data = list(sc,des3),method = "cali")
summary(fit)

fit = IPWM_build(data = list(sc,des3),method = "CLW")
summary(fit)

IPWM_estimate(fit,y='psa')
result = IPWM_estimate(fit,y='psa'); summary(result)

# with formula
fit = IPWM_build(data = list(sc,des1),method = "cali",
               p_formula = ~ agecat+marital+race+empstat+smoking+comorbidity+psa)
summary(fit)

# NA test
fit = IPWM_build(data = list(sc_na_test,des2),method = "ALP",
                 p_formula = ~ agecat+marital+race+empstat+smoking+comorbidity+height+weight,
                 na.action = na.exclude, sc_wname = "pseudo_wts")
head(fit$sc_updated)

head(sc_na_test)


# multi
fit_multi = IPWM_build(data = list(sc,des1,des2), sp_order = "given", verbose = T)
summary(fit_multi)
IPWM_estimate(fit_multi,'psa')

fit_multi = IPWM_build(
  data = list(sc,des1,des2),
  p_formula = list(~ psa, ~ agecat+marital+ empstat+ smoking+ comorbidity),
  precali = F,
  sp_order = 'size'
)



fit_multi = IPWM_build(data = list(sc,des1,des2),sp_order = 'given')
summary(fit_multi)

fit_multi = IPWM_build(data = list(sc,LFS=des1,CTADS=des2))
summary(fit_multi)

fit_multi = IPWM_build(
  data = list(sc,des1,des2),
  p_formula = list(~ agecat+marital, ~ empstat+ smoking+ comorbidity),
  precali = FALSE
)
summary(fit_multi)


fit_multi = IPWM_build(
  data = list(sc,des3,des2),
  p_formula = list(~ psa, ~ agecat+marital+race+
                    empstat+ smoking+ comorbidity+ height+weight),
  precali = TRUE,
  sp_order  = "size",
  na.action = na.omit,
  sc_wname = "wts"
)
summary(fit_multi)


# check survey.lonely.psu
library(survey)
?svyCprod
?surveyoptions
?svydesign
?svyrecvar


# lonely psu checking
getOption("survey.lonely.psu")
options(survey.lonely.psu = "fail")
options(survey.lonely.psu = "remove")
options(survey.lonely.psu = "adjust")
fit = IPWM_build(data = list(sc,des2),method = "alp")
summary(fit)
IPWM_estimate(fit,"psa")

# test zcol (binary variable works okay either as factor or numeric in BUILD stage)
sc_zcol = sc; sc_zcol$empstat = as.numeric(sc_zcol$empstat)-1; str(sc_zcol)
fit = IPWM_build(data = list(sc_zcol,des2),method = "alp")
head(fit$sc_updated)
fit = IPWM_build(data = list(sc,des2),method = "alp")
head(fit$sc_updated)
IPWM_estimate(fit,"psa",zcol = "empstat")


# test NR Error
sp2_test = sp2
sp2_test$BMI = sp2_test$weight/(sp2_test$height)^2
sp2_test$BMI = sp2_test$weight+(sp2_test$height)
head(sp2_test)
des2_test <- svydesign(
  ids = ~psu_sp2,
  strata = ~strata_sp2,
  weights = ~wts_sp2,
  data = sp2_test,
  nest = TRUE
)
fit = IPWM_build(data = list(sc,des2_test),method = "cali")



mean(sc$psa)
sum(sc$psa * sc$true_wts)/sum(sc$true_wts)




