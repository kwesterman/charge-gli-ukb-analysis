library(tidyverse)
library(data.table)


all_phenos_df <- fread("../data/raw/ukb47880.csv", nrows=1000, data.table=FALSE)

# Covariates -------------------------------------------------------------

basic_covar_df <- all_phenos_df %>%
  select(
    id = eid,
    age = `21022-0.0`,
    sex = `31-0.0`,  # Self-reported; 0=Female, 1=Male
    bmi = `21001-0.0`
  ) %>%
  mutate(
    age_sq = age ** 2
  )

genetic_covar_df <- all_phenos_df %>%
  select(
    id = eid,
    contains("22009-0")  # Genetic PC fields
  ) %>%
  rename_all(function(nm) gsub("22009-0.", "gPC", nm))

covar_df <- left_join(basic_covar_df, genetic_covar_df, by="id")

# Exposures --------------------------------------------------------------------

smoking_df <- all_phenos_df %>%
  select(id = eid)

alcohol_df <- all_phenos_df %>%
  select(id = eid)

sleep_df <- all_phenos_df %>%
  select(id = eid)

depression_df <- all_phenos_df %>%
  select(id = eid)

pa_df <- all_phenos_df %>%
  select(id = eid)

education_df <- all_phenos_df %>%
  select(id = eid)

exposure_df <- smoking_df %>%
  left_join(alcohol_df, by="id") %>%
  left_join(sleep_df, by="id") %>%
  left_join(depression_df, by="id") %>%
  left_join(pa_df, by="id") %>%
  left_join(education_df, by="id")

# Outcomes ---------------------------------------------------------------------

ukb_biomarker_fields <- c(
  alt = 30620, alb = 30600, alp = 30610, apoA = 30630, apoB = 30640,
  ast = 30650, hscrp = 30710, Ca = 30680, chol = 30690, creatinine = 30700,
  cysC = 30720, bilirubin_dir = 30660, ggt = 30730, glu = 30740, hba1c = 30750,
  hdl = 30760, igf1 = 30770, ldl = 30780, lipA = 30790, oestradiol = 30800,
  phos = 30810, rheum_factor = 30820, shbg = 30830, tes = 30850,
  bilirubin_tot = 30840, protein_tot = 30860, tg = 30870, urate = 30880,
  urea = 30670, vitD = 30890
)

biomarker_df <- all_phenos_df %>%
  select(
    id = eid,
    # all_of(setNames(paste0(ukb_biomarker_fields, "-0.0"), 
    #                 names(ukb_biomarker_fields)))
  )

outcome_df <- biomarker_df

# Final processing -------------------------------------------------------------

final_pheno_df <- covar_df %>%
  left_join(exposure_df, by="id") %>%
  left_join(outcome_df, by="id")

# Need to remove individuals who revoked consent

write_csv(final_pheno_df, "test_phenos.csv")
