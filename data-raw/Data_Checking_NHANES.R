library(dplyr)
library(skimr)
library(janitor)
library(DataExplorer)
library(naniar)
library(tidyverse)

sp_NHANES <- read.csv("NHANES_Harmonization_07302025.csv")


  sp_NHANES <- sp_NHANES %>%
  select(
    SDDSRVYR, d_agecat, d_marital, d_race, d_educat, d_empstat,
    d_bmicat2, d_smoking, d_arthritis, d_bronchitis, d_liver_comorb,d_osteoporosis,
    d_comorbidity,d_p_pros_cancer , d_psa_level0, WTINT10YR
  )

  sc  <- sp_NHANES %>% filter(SDDSRVYR == 3) %>% select(-SDDSRVYR)
  sc  <- sc %>% filter(!is.na(d_bmicat2))  # drop 3 NA values
  sc[ , 1:(ncol(sc) - 2)] <- lapply(sc[ , 1:(ncol(sc) - 2)], factor)
  sc <- sc %>% rename(psa = d_psa_level0, weight = WTINT10YR)
  str(sc)



  sp1 <- sp_NHANES %>% filter(SDDSRVYR != 3) %>% select(-d_psa_level0,-SDDSRVYR)
  sp1[ , 1:(ncol(sp1) - 1)] <- lapply(sp1[ , 1:(ncol(sp1) - 1)], factor)
  sp1 <- sp1 %>% rename(wts = WTINT10YR)
  str(sp1)

  dim(sc) ;summary(sc)
  dim(sp1);summary(sp1)




sp_NHIS <- read.csv("NHIS_Harmonization_07302025.csv")
sp_NHIS <- sp_NHIS %>% filter(SAMPWEIGHT != 0)
  sp2 <- sp_NHIS %>%
  select(
    d_agecat, d_marital, d_region, d_race, d_educat, d_empstat,
    d_smoking,d_bmicat2, d_arthritis, d_comorbidity,d_p_pros_cancer, SAMPWEIGHT
  )
  sp2 <- sp2 %>% rename(wts = SAMPWEIGHT)
  sp2[ , 1:(ncol(sp2) - 1)] <- lapply(sp2[ , 1:(ncol(sp2) - 1)], factor)
  str(sp2)

  dim(sp2);summary(sp2)

  summary(sp_NHIS)



# Initial inspection
dim(df)
head(df, 10)
tail(df, 10)
str(df,list.len=99999)         # Structure and data types
glimpse(df)     # dplyr-style overview of structure

# Summary statistics
summary(df)     # Basic summary stats
skim(df)        # Detailed summary: distribution, missing, unique values, etc.

# Missing values
# missing values per column
missing_per_col <- colSums(is.na(df))
print(missing_per_col/dim(df)[[1]])

# Visualize missing values
plot_missing(df)

# missing values per row
qplot(rowSums(is.na(df)), geom="histogram", bins=30,
      xlab="Missing values per row", ylab="Count")

# Data type of each column
sapply(df, class)

# Number of unique values per column
unique_vals <- sapply(df, function(x) length(unique(x)))
print(unique_vals)

# Actual unique values per column
unique_values_list <- lapply(df, unique)
print(unique_values_list)


# Categorical variable frequencies (need factor type)
df = sp
cat_df <- df %>% select(where(~ is.character(.) | is.factor(.)))
for (col in names(cat_df)) {
  cat("\nTop frequencies for", col, ":\n")
  print(sort(table(cat_df[[col]]), decreasing = TRUE)[1:10])
}


