## TO DO:
## - Must exclude individuals without genetics data prior to phenotyping (use existing of genetic PCs?)


# Preliminaries ---------------------------------------------------------

library(tidyverse)
library(data.table)

winsorize <- function(x, SDs=6) {
  bounds <- mean(x, na.rm=TRUE) + c(-1, 1) * SDs * sd(x, na.rm=TRUE)
  print(paste0(sum(x < bounds[1], na.rm=TRUE), " values winsorized at the lower bound."))
  print(paste0(sum(x > bounds[2], na.rm=TRUE), " values winsorized at the upper bound."))
  case_when(
    x < bounds[1] ~ bounds[1],
    x > bounds[2] ~ bounds[2],
    TRUE ~ x
  )
}

basic_phenos_df <- fread("../data/raw/ukb47880.csv", nrows=10000, data.table=FALSE)
blood_assays_df <- fread("../data/raw/ukb51002.csv", nrows=10000, data.table=FALSE)
urine_assays_df <- fread("../data/raw/ukb50582.csv", nrows=10000, data.table=FALSE)

# Covariates -------------------------------------------------------------

basic_covar_df <- basic_phenos_df %>%
  select(
    id = eid,
    age = `21022-0.0`,
    sex = `31-0.0`,  # Self-reported; 0=Female, 1=Male
    bmi = `21001-0.0`
  ) %>%
  mutate(
    age_sq = age ** 2,
    sex = 1 - sex  # Switch coding to 0=Male, 1=Female
  )

genetic_covar_df <- basic_phenos_df %>%
  select(
    id = eid,
    contains("22009-0")  # Genetic PC fields
  ) %>%
  rename_with(function(nm) gsub("22009-0.", "gPC", nm))

covar_df <- left_join(basic_covar_df, genetic_covar_df, by="id")

# Exposures --------------------------------------------------------------------

smoking_df <- basic_phenos_df %>%
  select(id = eid)

alcohol_df <- basic_phenos_df %>%
  select(id = eid)

sleep_df <- basic_phenos_df %>%
  select(
    id = eid,
    age = `21022-0.0`, 
    sex = `31-0.0`,  # Coding shouldn't matter for regression adjustment
    sleep_dur = `1160-0.0`
  ) %>%
  mutate(sleep_dur = ifelse(sleep_dur %in% c(-1, -3), NA, sleep_dur)) %>%
  filter(!is.na(sleep_dur))
sleep_df$sleep_dur_resid <- resid(lm(sleep_dur ~ sex * age, data=sleep_df))
sleep_df <- sleep_df %>%
  mutate(stst = as.integer(sleep_dur_resid <= quantile(sleep_dur_resid, 0.2)),
	 ltst = as.integer(sleep_dur_resid >= quantile(sleep_dur_resid, 0.8))) %>%
  select(id, sleep_dur_resid, stst, ltst)

depression_df <- basic_phenos_df %>%
  select(id = eid)

pa_df <- basic_phenos_df %>%
  select(id = eid)

education_df <- basic_phenos_df %>%
  select(id = eid)

exposure_df <- smoking_df %>%
  left_join(alcohol_df, by="id") %>%
  left_join(sleep_df, by="id") %>%
  left_join(depression_df, by="id") %>%
  left_join(pa_df, by="id") %>%
  left_join(education_df, by="id")

# Outcomes ---------------------------------------------------------------------

sr_meds <- basic_phenos_df %>%
  select(
    id = eid,
    contains("6177-0.")  # Self-reported cholesterol/BP/insulin medication use
  ) %>%
  mutate(
    bp_med = rowSums(.[, paste0("6177-0.", 0:2)] == 2, na.rm=TRUE) >= 1,
    chol_med = rowSums(.[, paste0("6177-0.", 0:2)] == 1, na.rm=TRUE) >= 1
  )

bp_df <- basic_phenos_df %>%
  select(
    id = eid,
    dbp1 = `4079-0.0`, dbp2 = `4079-0.1`,
    sbp1 = `4080-0.0`, sbp2 = `4080-0.1`
  ) %>%
  inner_join(sr_meds, by="id") %>%
  mutate(
    dbp = (dbp1 + dbp2) / 2,
    sbp = (sbp1 + sbp2) / 2,
    dbp = ifelse(bp_med, dbp + 10, dbp),
    sbp = ifelse(bp_med, sbp + 15, sbp),
    pp = sbp - dbp
  ) %>%
  mutate(across(c(dbp, sbp, pp), winsorize)) %>%
  select(id, dbp, sbp, pp, bp_med)


ukb_biomarker_fields <- c(
  chol = 30690,
  glu = 30740, hba1c = 30750,
  hdl = 30760, ldl = 30780,
  tg = 30870
)

biomarker_df <- blood_assays_df %>%
  select(
    id = eid,
    all_of(setNames(paste0(ukb_biomarker_fields, "-0.0"),
                   names(ukb_biomarker_fields)))
  )

lipids_df <- biomarker_df %>%
  select(id, hdl, tg, ldl) %>%
  # Need fasting adjustment here
  inner_join(sr_meds, by="id") %>%
  mutate(ldl = ifelse(chol_med, ldl / 0.7, ldl)) %>%
  mutate(across(c(hdl, tg), log)) %>%  # Confirm that log-transform should happen AFTER meds adjustment
  mutate(across(c(hdl, tg, ldl), winsorize)) %>%
  select(id, hdl, tg, ldl)
  
outcome_df <- full_join(
  bp_df, 
  lipids_df,
  by="id"
)

# Final processing -------------------------------------------------------------

final_pheno_df <- covar_df %>%
  left_join(exposure_df, by="id") %>%
  left_join(outcome_df, by="id")

# Need to remove individuals who revoked consent

write_csv(final_pheno_df, "test_phenos.csv")
