## code to prepare `DATASET` dataset goes here

library(tidyverse)

# update .Rdata
save(list = ls(), file = "C:/Users/Kevin/OneDrive - University of Toronto/Summer_Job/Github_IWH/Package6 - BackUp - D matrix and Seperation/IPWMethod/.RData")


# PREPARE sp_NHANES
sp_NHANES <- read.csv("NHANES_Harmonization_07302025.csv")
summary(sp_NHANES)
str(sp_NHANES)


sum(sp_NHANES$WHD010 == 9999, na.rm = TRUE)
sum(sp_NHANES$WHD020 == 99999, na.rm = TRUE)
sum(sp_NHANES$WHD020 == 9999, na.rm = TRUE)

sp_NHANES$WHD010[sp_NHANES$WHD010 == 9999] <- NA
sp_NHANES$WHD020[sp_NHANES$WHD020 == 99999] <- NA
sp_NHANES$WHD020[sp_NHANES$WHD020 == 9999] <- NA

sp_NHANES <- sp_NHANES %>%
  select(
    SDDSRVYR, d_agecat, d_marital, d_race, d_educat, d_empstat,
    d_smoking, d_comorbidity , d_psa_level0, WHD010, WHD020, WTINT10YR,
    SDMVSTRA, SDMVPSU
  )
sp_NHANES <- sp_NHANES %>% rename( agecat = d_agecat, marital = d_marital, race = d_race, educat = d_educat,
                                   empstat = d_empstat, smoking = d_smoking, comorbidity = d_comorbidity,
wave = SDDSRVYR , psa = d_psa_level0, height = WHD010, weight = WHD020, true_wts =WTINT10YR,
psu  = SDMVPSU, strata= SDMVSTRA
)


sp_NHANES <- sp_NHANES %>%
  mutate(
    height = height * 0.0254,
    weight = weight * 0.453592,
    BMI = weight / (height^2)
  )
n <- ncol(sp_NHANES)
sp_NHANES <- sp_NHANES[, c(1:(n-2), n, n-1)]


summary(sp_NHANES)



sc_raw  <- sp_NHANES %>% filter(wave %in% c(2, 3)) %>% select(-wave,-psu,- strata)

sc_raw[ , 1:(ncol(sc_raw) - 5)] <- lapply(sc_raw[ , 1:(ncol(sc_raw) - 5)], factor)
summary(sc_raw)

sc_clean  <- na.omit(sc_raw)
summary(sc_clean)
sc = sc_clean



sp1_raw <- sp_NHANES %>% filter(wave %in% c(4, 5)) %>% select(-wave,-height,-weight,-BMI)  %>% rename ( wts_sp1 = true_wts )
summary(sp1_raw)
sp1_raw[ , 1:(ncol(sp1) - 4)] <- lapply(sp1_raw[ , 1:(ncol(sp1) - 4)], factor)
sp1_clean = na.omit(sp1_raw)
sp1 = sp1_clean
str(sp1)




# sp2

sp_NHIS <- read.csv("NHIS_Harmonization_07302025.csv") %>% filter(SAMPWEIGHT != 0)


sp2 <- sp_NHIS %>%
  select(d_agecat, d_marital, d_race, d_empstat,
         d_bmicat2, d_smoking, d_comorbidity , HEIGHT,WEIGHT ,SAMPWEIGHT,
         STRATA,PSU)

sp2 <- sp2 %>% rename(  agecat = d_agecat, marital = d_marital, race = d_race, empstat = d_empstat,
                        bmicat = d_bmicat2, smoking = d_smoking, comorbidity = d_comorbidity,
  height=HEIGHT , weight=WEIGHT , wts_sp2 = SAMPWEIGHT, strata_sp2=STRATA, psu_sp2=PSU)





sp2$height[sp2$height == 99] <- NA
sp2$weight[sp2$weight == 999] <- NA




sp2 <- sp2 %>%
  mutate(
    height = height * 0.0254,
    weight = weight * 0.453592,
    BMI = weight / (height^2)
  )
dim(sp2);summary(sp2)
n= ncol(sp2)
sp2 <- sp2[, c(1:(n-2), n, n-1)]



dim(sp2);summary(sp2)
sp2[ , 1:(ncol(sp2) - 6)] <- lapply(sp2[ , 1:(ncol(sp2) - 6)], factor)
sp2 = na.omit(sp2)

summary(sp2)

sc_full = sc_raw; sc_cc = sc_clean; sp1_full = sp1_raw; sp1_cc = sp1_clean; sp2_full = sp2_raw; sp2_cc = sp2_clean
sc = sc_cc; sp1 = sp1_cc; sp2 = sp2_cc




sc = sc_cc
ref_survey_1 = sp1_cc
ref_survey_2 = sp2_cc
rm(sp1,sp2)

dim(sc);str(sc)
dim(ref_survey_1);str(ref_survey_1)
dim(ref_survey_2);str(ref_survey_2)


# create .rda in data file
sc_with_na = sc_na_test
ref_survey_1_psu = sp1
ref_survey_2 = sp2
ref_survey_1_bootstrap = sp1_bootstrap

usethis::use_data(sc, overwrite = TRUE)
usethis::use_data(sc_with_na, overwrite = TRUE)

usethis::use_data(ref_survey_1_psu, overwrite = TRUE)
usethis::use_data(ref_survey_1_bootstrap, overwrite = TRUE)
usethis::use_data(ref_survey_2, overwrite = TRUE)



# write documentation for data
usethis::use_r("sc")
usethis::use_r("ref_survey_1")
usethis::use_r("sef_survey_2")

