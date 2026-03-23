# -------- Setting Up ------------
options(scipen = 999) #avoid sci expression
set.seed(123456)

N=500000

v1 = rbinom(N, size=1, prob=0.5)                # Bernoulli(0.5)
v2 = runif(N, min=0, max=2)                     # Uniform(0, 2)
v3 = rexp(N, rate = 1)                          # Exponential(1)
v4 = rchisq(N, df = 4)                          # Chi-square with 4 degrees of freedom

# Covariates
x1 = v1
x2 = v2 #+ 0.3 * x1
x3 = v3 #+ 0.2 * (x1 + x2)
x4 = v4 #+ 0.1 *(x1 + x2 + x3)

mu =  - x1 -  x2 + x3 + x4
y = mu + rnorm(N)
mu_true = mean(y)


np=12500

a <- x3 + 0.03*y
rng <- range(a)
if (diff(rng) == 0) stop("x3 + 0.03*y is constant; no solution (ratio is always 1).")
amin <- rng[1]; amax <- rng[2]

cnst_sp <- (amax - 20*amin) / 19

# sanity check
q <- cnst_sp + a
(max(q) / min(q))  # should be 20

pi_p1=np * q / sum(q)  #inclusion probability for fp
max(pi_p1)/min(pi_p1)   # 20

# pi_p are all valid
any(is.na(pi_p1))        #False
all(pi_p1>0 & pi_p1<1)    #True
di1=1/pi_p1




np=25000

a <- x3 + 0.03*y
rng <- range(a)
if (diff(rng) == 0) stop("x3 + 0.03*y is constant; no solution (ratio is always 1).")
amin <- rng[1]; amax <- rng[2]

cnst_sp <- (amax - 20*amin) / 19

# sanity check
q <- cnst_sp + a
(max(q) / min(q))  # should be 20

pi_p2=np * q / sum(q)  #inclusion probability for fp
max(pi_p2)/min(pi_p2)   # 20

# pi_p2 are all valid
any(is.na(pi_p2))        #False
all(pi_p2>0 & pi_p2<1)    #True
di2=1/pi_p2



fp = data.frame( x1 = x1, x2 = x2, x3 = x3, x4 = x4 , y = y, wt_sp1=di1, wt_sp2=di2 )


nc = 2500
eta <-  0.18 * x1 + 0.18 * x2 - 0.27 * x3 - 0.27 * x4
exp_beta0 = nc/(sum(exp(eta)))
beta0=log(exp_beta0)

# verify       sum(exp(eta+beta0))    #30000

pi_c  =  exp(eta+beta0)
all(pi_c>0 & pi_c<1)
fp$wt_sc=1/pi_c



# --- numeric datasets --- #
sc   = fp[rbinom(N, size = 1, pi_c) == 1,c('x1','x2','x3','x4','y')]
sp   = fp[rbinom(N, 1, pi_p1) == 1,c('x1','x2','x3','x4', 'wt_sp1')]
sp1  = fp[rbinom(N, 1, pi_p1) == 1,c('x1','x2', 'wt_sp1')]
sp2  = fp[rbinom(N, 1, pi_p2) == 1,c('x3','x4', 'wt_sp2')]

sp3  = fp[rbinom(N, 1, pi_p1) == 1,c('x1','x2','x3', 'wt_sp1')]
sp4  = fp[rbinom(N, 1, pi_p2) == 1,c('x3','x4','x2', 'wt_sp2')]


# --- categorical datasets --- #

# For sc
x3_factor <- factor(cut(sc$x3,
                        breaks = quantile(sc$x3, probs = seq(0, 1, by = 0.25), na.rm = TRUE),
                        include.lowest = TRUE,
                        labels = FALSE))

x4_factor <- factor(cut(sc$x4,
                        breaks = quantile(sc$x4, probs = seq(0, 1, by = 0.25), na.rm = TRUE),
                        include.lowest = TRUE,
                        labels = FALSE))

sc_cate <- data.frame(
  x1 = sc$x1,
  x2 = sc$x2,
  x3 = x3_factor,
  x4 = x4_factor,
  y  = sc$y
)

# For sp
x3_factor <- factor(cut(sp$x3,
                        breaks = quantile(sp$x3, probs = seq(0, 1, by = 0.25), na.rm = TRUE),
                        include.lowest = TRUE,
                        labels = FALSE))

x4_factor <- factor(cut(sp$x4,
                        breaks = quantile(sp$x4, probs = seq(0, 1, by = 0.25), na.rm = TRUE),
                        include.lowest = TRUE,
                        labels = FALSE))

sp_cate <- data.frame(
  x1 = sp$x1,
  x2 = sp$x2,
  x3 = x3_factor,
  x4 = x4_factor,
  wt_sp = sp$wt_sp1
)

PRECALI(c('y','wt_sp'),list(sc_cate,sp_cate),dup_vars = c("x3","x4"))

library(survey)





# always call raking "Calibration"
# Pseudo-weighted (multi) ALP CLW Calibration

# --------------- IPWM Showcase --------------------

# one-ref case
head(sc);head(sp)

IPWM(sc=sc,sp=sp,y='y',weight = 'wt_sp1',method = "cali")

outcome = IPWM(sc=sc,sp=sp,y='y',weight = 'wt_sp1',method = "calibration")
summary(outcome)


# wrong use of p_formula
outcome = IPWM(sc=sc,sp=sp,y='y',weight = 'wt_sp1',method = "AMP", p_formula = c('x1','x3'))
summary(outcome)

# wrong y
outcome = IPWM(sc=sc,sp=sp,y='y',weight = 'wt_sp1',method = "ALP", p_formula = ~ x1 + x5)
summary(outcome)

# wrong weight column
outcome = IPWM(sc=sc,sp=sp,y='y',weight = 'wt_sp2',method = "ALP", p_formula = ~ x1 + x3)
summary(outcome)


outcome = IPWM(sc=sc,sp=sp,y='y',weight = 'wt_sp1',method = "raking", p_formula = ~ x1 + x2 + x3 + x4+ I(x3^2))
summary(outcome)

outcome = IPWM(sc=sc,sp=sp,y='y',weight = 'wt_sp1',method = "ALP", p_formula = ~ x1 + x2 + x3 * x4)
summary(outcome)


# categorical
head(sc_cate);head(sp_cate)

# Basic additive models
outcome = IPWM(sc=sc_cate,sp=sp_cate,y='y',weight = 'wt_sp',method = "ALP", p_formula = ~ x1 + x2 + x3 + x4)
summary(outcome)

outcome = lm(y ~ x1 + x2 + x3 + x4,data = sc_cate)

# Interaction between two factors
outcome = IPWM(sc=sc_cate,sp=sp_cate,y='y',weight = 'wt_sp',method = "ALP", p_formula = ~ x1 * x2 + x3 * x4)
summary(outcome)

# Factor + numeric interaction
outcome = IPWM(sc=sc_cate,sp=sp_cate,y='y',weight = 'wt_sp',method = "ALP", p_formula = ~ x1 * x3 + x2 + x4)
summary(outcome)

outcome = IPWM(sc=sc_cate,sp=sp_cate,y='y',weight = 'wt_sp',method = "ALP", p_formula = NULL)
summary(outcome)

# Full interaction between all predictors
summary(outcome)




# multi-ref case
IPWM(sc=sc,sp=list(sp1,sp2),y='y', weight = c('wt_sp1','wt_sp2'), cali = T)
IPWM(sc=sc,sp=list(sp3,sp4),y='y', weight = c('wt_sp1','wt_sp2'), cali = T)
IPWM(sc=sc,sp=list(sp3,sp4),y='y', weight = c('wt_sp1','wt_sp2'), cali = T,method = "ALP")




result = IPWM(sc=sc,sp=list(sp2,sp1),y='y', weight = c('wt_sp2','wt_sp1'), cali = T)
summary(result)

head(sp3);head(sp4)
result = IPWM(sc=sc,sp=list(sp3,sp4),y='y', weight = c('wt_sp1','wt_sp2'), cali = T)
summary(result)

result = IPWM(sc=sc,sp=list(CTADS=sp3,LFS=sp4),y='y', weight = c('wt_sp1','wt_sp2'), cali = T)
summary(result)

CTADS = sp1; LFS=sp2
sp_list = list(CTADS = sp1, LFS=sp2);names(sp_list)

result = IPWM(sc=sc,sp = list(CTADS = sp1, LFS = sp2),y='y', weight = c('wt_sp1','wt_sp2'), cali = T, p_formula = list( ~ x1 , ~ x3 + x4))
summary(result)

result = IPWM(sc=sc,sp=list(sp3,sp4),y='y', weight = c('wt_sp1','wt_sp2'), cali = T, p_formula = list(~ x1*x2  , ~ x2*x3 ))
summary(result)


result = IPWM(sc=sc,sp=list(CTADS = sp3,LFS =sp4),y='y', weight = c('wt_sp1','wt_sp2'), cali = T, p_formula = list(~ x1 + x2  , ~x3 + x4))
summary(result)


result = IPWM(sc=sc,sp=list(CTADS = sp1, LFS=sp2),y='y', weight = c('wt_sp1','wt_sp2'), cali = T, p_formula = list(~ x1, ~x3))
summary(result)

result = IPWM(sc=sc,sp=list(sp1,sp2),y='y', weight = c('wt_sp1','wt_sp2'), cali = T, p_formula = list(~ x1, ~ I(x3^2)))
summary(result)

result = IPWM(sc=sc,sp=list(sp1,sp2),y='y', weight = c('wt_sp1','wt_sp2'), cali = T, p_formula = list(~ x1 + x2, ~ x3 * x4 ))
summary(result)

result = IPWM(sc=sc,sp=list(sp1,sp2),y='y', weight = c('wt_sp1','wt_sp2'), cali = T, p_formula = list(~ x1 + x2, ~ x3 * x4 ))
summary(result)

result <- IPWM(
  sc = sc,
  sp = list(CTADS = sp1, LFS = sp2),
  y = 'y',
  weight = c('wt_sp1', 'wt_sp2'),
  cali = TRUE,
  p_formula = list(~ x1, ~ x3 + x4)
)
head(sp2)





# real data example




data(sc)


summary(sc)
sc <- sc %>% filter(!is.na(psa))


data("ref_survey_1")
summary(ref_survey_1)
ref_survey_1 <- ref_survey_1 %>% filter(!is.na(d_bmicat2))

result = IPWM(sc=sc, sp=ref_survey_1,  y='psa', weight = 'wts', method = "ALP")
result = IPWM(sc=sc, sp=ref_survey_1,  y='psa', weight = 'wts_sp1', p_formula = ~d_marital+d_agecat, method = "ALP")
result = IPWM(sc=sc, sp=ref_survey_1,  y='psa', weight = 'wts', p_formula = ~d_bmicat2+d_agecat+d_marital+d_race+d_educat+d_empstat+
                d_smoking+d_comorbidity, method = "ALP")



result = IPWM(sc=sc, sp=sp1,  y='psa', weight = 'wts_sp1', method = "ALP")

result = IPWM(sc=sc, sp=sp1,  y='psa', weight = 'wts', p_formula = ~d_bmicat2+d_agecat+d_marital+d_race+d_educat+d_empstat+
                d_smoking+d_comorbidity, method = "ALP")


result = IPWM(sc=sc, sp=list(ref1=sp1,ref2=sp2), y='psa', weight = c('wts','wts'), p_formula = list(~d_bmicat2+d_comorbidity,~d_smoking+d_empstat))


result = IPWM(sc=sc_clean, sp=list(sp1_clean,sp2_clean), y='psa', weight = c('wts_sp1','wts_sp2'), p_formula = list( ~~d_agecat + d_marital + d_race  + d_empstat   +    d_smoking + d_comorbidity , ~ height + weight ))
summary(result)



summary(result)

result = IPWM(sc=sc, sp=ref_survey_1, y='psa', weight = 'wts',method = "CLW")
summary(result)

result = IPWM(sc=sc, sp=ref_survey_1, y='psa', weight = 'wts',method = "calibration")
summary(result)


result = IPWM(sc=sc, sp=ref_survey_1, y='psa', weight = 'wts', p_formula = ~ d_agecat + d_marital + d_race +
              d_educat + d_empstat + d_bmicat2+ d_smoking  +  d_comorbidity, method = "calibration")
summary(result)


# with domain variable
dummy_vars <- model.matrix(~ d_bmicat2 - 1, data = sc)
head(dummy_vars)

# Combine sc with the dummy vars
sc <- cbind(sc, dummy_vars)
head(sc)

result = IPWM(sc=sc, sp=ref_survey_1, y='psa', weight = 'wts' ,zcol = "d_bmicat21" ,method = "calibration")
summary(result)

result = IPWM(sc=sc, sp=ref_survey_1, y='psa', weight = 'wts' ,zcol = "d_bmicat22" ,method = "calibration")
summary(result)

result = IPWM(sc=sc, sp=ref_survey_1, y='psa', weight = 'wts' ,zcol = "d_bmicat23" ,method = "calibration")
summary(result)

result = IPWM(sc=sc, sp=ref_survey_1, y='psa', weight = 'wts' ,zcol = "d_bmicat24" ,method = "calibration")
summary(result)


# true weighted psa mean
head(sc)
result <- sc %>%
  group_by(d_bmicat2) %>%
  summarise(
    weighted_mean_psa = sum(psa * weight) / sum(weight),
    n = n()
  )
result  # 2.04 1.67 1.56 2.04

# true weighted psa mean
head(sc)
result <- sc %>%
  summarise(
    weighted_mean_psa = sum(psa * weight) / sum(weight),
    n = n()
  )
result  # 2.04 1.67 1.56 2.04


# two sample case
data("ref_survey_2")
summary(ref_survey_2)

help(IPWM)

result = IPWM(sc=sc_cc, sp=list(sp1_cc,sp2_cc), y='psa', weight = c('wts_sp1','wts_sp2'),  cali = TRUE,
              p_formula = list( ~ psa ,
                                ~ agecat + marital + race + empstat + smoking + comorbidity +  BMI))
summary(result)


summary(result)
result = IPWM(sc=sc_full, sp=list(sp1_full,sp2_full), y='psa', weight = c('wts_sp1','wts_sp2'),  cali = TRUE,
              p_formula = list( ~ psa ,
                                ~ agecat + marital + race + empstat + smoking + comorbidity +  BMI))
summary(result)
summary(sp2_cc)


# --- test cumulative calibration --- #


# --- 0) build SP1/SP2/SP3 as defined ---
## SP1: first 500，agecat, marital, wts1
SP1 <- sp1_cc[1:400, c("agecat", "marital", "wts_sp1")]
names(SP1)[names(SP1) == "wts_sp1"] <- "wts1"

## SP2: 501–1000，agecat, race, wts2
SP2 <- sp1_cc[401:1000, c("agecat", "race", "wts_sp1")]
names(SP2)[names(SP2) == "wts_sp1"] <- "wts2"

## SP3: 1001–1240，agecat, marital, race, educat, wts3
SP3 <- sp1_cc[1001:1240, c("agecat", "marital", "race", "educat", "wts_sp1")]
names(SP3)[names(SP3) == "wts_sp1"] <- "wts3"


# --- 1) expand each SP to model matrix (no intercept), then bind weights back ---
expand_to_df <- function(df, wname) {
  vars <- setdiff(names(df), wname)

  # important: keep factors as factors (if they are already), model.matrix will do dummies
  X <- stats::model.matrix(~ ., data = df[, vars, drop = FALSE])
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  # make it a data.frame + add weight column
  out <- as.data.frame(X)
  out[[wname]] <- df[[wname]]
  out
}
SP1_mm <- expand_to_df(SP1, "wts1")
SP2_mm <- expand_to_df(SP2, "wts2")
SP3_mm <- expand_to_df(SP3, "wts3")

# --- 2) prepare sp_raw / sp_new / weight list ---
sp_raw <- list(SP1 = SP1_mm, SP2 = SP2_mm, SP3 = SP3_mm)

# initial sp_new can be the same (since weights are the starting weights)
sp_new <- sp_raw

weight <- list("wts1", "wts2", "wts3")

# --- 3) run cumulative pre-calibration ---
fit <- precal_cumulative(sp_raw = sp_raw, sp_new = sp_new, weight = weight, verbose = TRUE)

# --- 4) get weighted totals in model matrix
weighted_totals_mm <- function(df, wname) {
  X_cols <- setdiff(names(df), wname)
  w <- df[[wname]]

  c(
    "(Intercept)" = sum(w),
    colSums(as.matrix(df[, X_cols, drop = FALSE]) * w)
  )
}


SP1_cal <- fit$sp_new$SP1
SP2_cal <- fit$sp_new$SP2
SP3_cal <- fit$sp_new$SP3

tot_SP1 <- weighted_totals_mm(SP1_cal, "wts1")
tot_SP1
tot_SP2 <- weighted_totals_mm(SP2_cal, "wts2")
tot_SP2
tot_SP3 <- weighted_totals_mm(SP3_cal, "wts3")
tot_SP3




result = IPWM(sc=sc, sp=list(SP1,SP2,SP3), y='psa',
              weight = c('wts1','wts2','wts3'),  precali = TRUE, sp_order = "size")

result = IPWM(sc=sc, sp=list(SP1,SP2), y='psa',
              weight = c('wts1','wts2'),  precali = TRUE, sp_order = 'given')

result = IPWM(sc=sc, sp=list(SP1,SP2,SP3), y='psa',
              weight = c('wts1','wts2','wts3'),  precali = FALSE, sp_order = "size")

result = IPWM(sc=sc, sp=list(SP1,SP2,SP3), y='psa',
              weight = c('wts1','wts2','wts3'),  precali = FALSE, sp_order = "given")

summary(result)














# --------------- IPWM Showcase with design ----------------

save(list = ls(), file = "C:/Users/Kevin/OneDrive - University of Toronto/Summer_Job/Github_IWH/Package6 - BackUp - D matrix and Seperation/IPWMethod/.RData")

head(sc);head(sp1);head(sp2)

attr(sc, "na.action") <- NULL
attr(sp1, "na.action") <- NULL

outcome = IPWM(sc=sc,sp=sp1,y='psa',weight = 'wts_sp1',method = "calibration",
               p_formula = ~ agecat+marital+race+educat+empstat+smoking+comorbidity,
               sp_design = list(type='poisson'))
summary(outcome)

outcome = IPWM(sc=sc,sp=sp1,y='psa',weight = 'wts_sp1',method = "calibration", sp_design = list(type='pps'))
summary(outcome)

outcome = IPWM(sc=sc,sp=sp1,y='psa',weight = 'wts_sp1',method = "calibration",
                       sp_design = list(
                         type   = "psu",
                         ids    = "psu",      # PSU column name in sp
                         strata = "strata",   # strata column name in sp
                         fpc    = NULL,       # optional: column name for fpc, or NULL
                         nest   = TRUE      ))
summary(outcome)

outcome = IPWM(sc=sc,sp=sp1_bt,y='psa',weight = 'wts_sp1',method = "calibration",
             sp_design = list(
             type       = "replicate",
             weights    = "wts_sp1",                 # main weight column name in sp
             repweights = paste0("bw", 1:500), # replicate weight column names in sp
             rep_type   = "bootstrap"          # other options to be added
             ))
summary(outcome)


## SP2
head(sc);head(sp2)
attr(sp2, "na.action") <- NULL
outcome = IPWM(sc=sc,sp=sp2,y='psa',weight = 'wts_sp2',method = "calibration",
               p_formula = ~ agecat+marital+race+empstat+smoking+comorbidity+BMI,
               sp_design = list(
                 type   = "psu",
                 ids    = "psu_sp2",      # PSU column name in sp
                 strata = "strata_sp2",   # strata column name in sp
                 fpc    = NULL,       # optional: column name for fpc, or NULL
                 nest   = TRUE      ))
summary(outcome)




head(ref_survey_2)
sp2 = ref_survey_2

# multi case
outcome = IPWM(sc=sc, sp= list(sp1,sp2), y='psa',
               weight = c('wts_sp1',"wts_sp2") ,  precali = T, sp_order = 'size',
               p_formula = list(~ psa,
                                ~ agecat+marital+race+empstat+smoking+comorbidity+BMI),
              sp_design = list( list(
                                type   = "psu",
                                ids    = "psu",      # PSU column name in sp
                                strata = "strata",   # strata column name in sp
                                nest   = TRUE      ),
                                list(
                                  type   = "psu",
                                  ids    = "psu_sp2",      # PSU column name in sp
                                  strata = "strata_sp2",   # strata column name in sp
                                  nest   = TRUE      )))

IPWM(sc=sc,sp=list(sp1,sp2),y="psa",weight=c("wts_sp1","wts_sp2"),
     precali=TRUE,sp_order="size",
     p_formula=list(~psa,~agecat+marital+race+empstat+smoking+comorbidity+BMI),
     sp_design = list(list(type='poisson'),list(type='poisson')))

head(sc)
sum(sp1$psa * sp1$wts_sp1) /sum(sp1$wts_sp1)
summary(outcome)




