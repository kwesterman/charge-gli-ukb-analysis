# Parse command-line argument(s) -----------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
tag <- ifelse(is.na(args[1]), "test", args[1]) # Usually should be the date in YYYYMMDD format

# Input filenames --------------------------------------------------------------

basic_phenos_file <- "../data/raw/ukb47880.csv"
blood_assays_file <- "../data/raw/ukb51002.csv"
urine_assays_file <- "../data/raw/ukb50582.csv"
pack_years_file <- "../data/raw/ukb668908.csv"

med_atc_file <- "../data/raw/wu_et_al/41467_2019_9572_MOESM3_ESM_tw.txt"
med_atc_file <- "../data/raw/41467_2019_9572_MOESM3_ESM_tw.txt"

pan_ancestry_file <- "../data/raw/Files\ for\ retman/all_pops_non_eur_pruned_within_pop_pc_covs.tsv"
pa_bridge_file <- "../data/raw/ukb8343bridge31063.txt"
sample_exclusion_file <- "../data/raw/w8343_20220222.csv"

# Preliminaries ----------------------------------------------------------------

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

pa_bridge_df <- fread(pa_bridge_file, col.names=c("id", "s"), data.table=FALSE)
pan_ancestry_df <- fread(pan_ancestry_file, data.table=FALSE) %>%
  inner_join(pa_bridge_df, by="s")

sample_exclusions <- readLines(sample_exclusion_file)  # Participants who revoked consent

valid_ids <- setdiff(pan_ancestry_df$id, sample_exclusions)  # Participants with genetics but not having revoked consent

pack_years_df <- fread(pack_years_file, nrows=10000, data.table=FALSE) # Load pack years first to merge into basic pheno data
basic_phenos_df <- fread(basic_phenos_file, nrows=Inf, data.table=FALSE) %>%
  left_join(pack_years_df, by="eid") %>%
  filter(eid %in% valid_ids)
blood_assays_df <- fread(blood_assays_file, nrows=10000, data.table=FALSE) %>%
  filter(eid %in% valid_ids)
urine_assays_df <- fread(urine_assays_file, nrows=10000, data.table=FALSE) %>%
  filter(eid %in% valid_ids)
med_atc_df <- fread(med_atc_file, data.table=FALSE) 

# Covariates -------------------------------------------------------------------

basic_covar_df <- basic_phenos_df %>%
  select(
    id = eid,
    age = `21022-0.0`,
    sex = `31-0.0`,  # Self-reported; 0=Female, 1=Male
    bmi = `21001-0.0`,
    ethnicity = `21000-0.0`,
    assessment_center = `54-0.0`,
    fasting_time = `74-0.0`
  ) %>%
  mutate(
    age_sq = age ** 2,
    sex = 1 - sex,  # Switch coding to 0=Male, 1=Female
    ethnicity = ifelse(ethnicity %in% c(-1, -3), NA, ethnicity)
  )

genetic_covar_df <- pan_ancestry_df %>%
  select(id, pop, related, contains("PC")) %>%
  rename_with(function(nm) sub("^PC", "gPC", nm))

covar_df <- left_join(basic_covar_df, genetic_covar_df, by="id")

# Exposures --------------------------------------------------------------------

# Not currently calculating qCPD or dPY25 since they were removed from the analysis plans
# Analysis plan states to exclude current smokers with CPD or PY = 0 (they are set to NA below)
smoking_df <- basic_phenos_df %>%
  select(id = eid,
         age = `21022-0.0`,
         smoker_status = `20116-0.0`,
         cigs_per_day = `3456-0.0`,
         pack_years = `20161-0.0`,
         age_start_smoking = `3436-0.0`,
         age_stop_smoking = `2897-0.0`, # age stopped smoking cigarettes
         ever_stop_6mo = `2907-0.0`,
         ever_try_stop = `3486-0.0`
        ) %>%
  mutate(cursmk = ifelse(is.na(smoker_status) | smoker_status == -3, NA,
                         ifelse(smoker_status == 2, 1, 0)),
         cpd = as.numeric(ifelse((cigs_per_day %in% c(-1, -3, -10, 0)) | is.na(cursmk) | cursmk==0, NA, cigs_per_day)), # -10 = < 1 per day
         py = ifelse((pack_years %in% c(-1, -3, -10, 0)) | is.na(cursmk) | cursmk==0, NA, pack_years)
        ) %>%
  select(id, cursmk, cpd, py)

alcohol_df <- basic_phenos_df %>%
  select(id = eid,
         sex = `31-0.0`,  # Self-reported; 0=Female, 1=Male
         drinker_status = `20117-0.0`,
         n_wk_beer = `1588-0.0`,
         n_wk_wine_red = `1568-0.0`,
         n_wk_wine_wht_chm = `1578-0.0`,
         n_wk_wine_fort = `1608-0.0`,
         n_wk_spirits = `1598-0.0`,
         n_wk_other = `5364-0.0`
        ) %>%
  mutate(sex = 1 - sex, # Switch coding to 0=Male, 1=Female
         n_wk_beer = ifelse(is.na(n_wk_beer) | n_wk_beer %in% c(-1, -3), 0, n_wk_beer),
         n_wk_wine_red = ifelse(is.na(n_wk_wine_red) | n_wk_wine_red %in% c(-1, -3), 0, n_wk_wine_red),
         n_wk_wine_wht_chm = ifelse(is.na(n_wk_wine_wht_chm) | n_wk_wine_wht_chm %in% c(-1, -3), 0, n_wk_wine_wht_chm),
         n_wk_wine_fort = ifelse(is.na(n_wk_wine_fort) | n_wk_wine_fort %in% c(-1, -3), 0, n_wk_wine_fort),
         n_wk_spirits = ifelse(is.na(n_wk_spirits) | n_wk_spirits %in% c(-1, -3), 0, n_wk_spirits),
         n_wk_other = ifelse(is.na(n_wk_other) | n_wk_other %in% c(-1, -3), 0, n_wk_other),
         curdrink = ifelse(is.na(drinker_status) | drinker_status == -3, NA,
                           ifelse(drinker_status == 2, 1, 0)),
         never = ifelse(!is.na(drinker_status) & drinker_status == 0, 1, 0),
         total_alc_uk = (n_wk_beer * 16 +
                         n_wk_wine_red *16.8 +
                         n_wk_wine_wht_chm * 16.8 +
                         n_wk_wine_fort * 14.08 +
                         n_wk_spirits * 8 +
                         n_wk_other * 12),
         drinks_per_week = total_alc_uk / 14 # Weekly total US standard measure = Weekly total UK g of alcohol / 14 g per measure
        ) %>%
  mutate(drinks_per_week = ifelse(is.na(curdrink) | curdrink == 0, NA, drinks_per_week),
         light = ifelse(!is.na(curdrink) & curdrink == 1 &
                        ((sex == 0 & drinks_per_week <= 14) |
                         (sex == 1 & drinks_per_week <= 7)), 1, 0),
         heavy = ifelse(!is.na(curdrink) & curdrink == 1 &
                        ((sex == 0 & drinks_per_week > 14) |
                         (sex == 1 & drinks_per_week > 7)), 1, 0),
         very_heavy = ifelse(!is.na(curdrink) & curdrink == 1 &
                             ((sex == 0 & drinks_per_week > mean(drinks_per_week[sex==0],na.rm=TRUE) + 6*sd(drinks_per_week[sex==0],na.rm=TRUE)) |
                              (sex == 1 & drinks_per_week > mean(drinks_per_week[sex==1],na.rm=TRUE) + 6*sd(drinks_per_week[sex==1],na.rm=TRUE))), 1, 0),
         heavy_vs_light = ifelse(is.na(curdrink) | curdrink==0 | very_heavy, NA, ifelse(heavy, 1, 0)),
         light_vs_never = ifelse(!never & !light, NA, ifelse(light, 1, 0)),
         heavy_vs_never = ifelse(!never & !heavy, NA, ifelse(heavy, 1, 0))) %>%
  select(id, curdrink, heavy_vs_light, light_vs_never, heavy_vs_never)

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

ls_atc_codes = strsplit(med_atc_df$Medication_ATC_code," |",fixed=T)
is_bp_med = unlist(lapply(ls_atc_codes, function(x) any(grepl("^C02|^C03|^C07|^C08|^C09",x))))
is_chol_med = unlist(lapply(ls_atc_codes, function(x) any(grepl("^C10",x))))

meds_df <- basic_phenos_df %>%
  select(
    id = eid,
    contains("6177-0."),  # Self-reported cholesterol/BP/insulin medication use
	contains("20003-0.")
  ) %>%
  mutate(
    bp_med_sr = rowSums(.[, paste0("6177-0.", 0:2)] == 2, na.rm=TRUE) >= 1,
    chol_med_sr = rowSums(.[, paste0("6177-0.", 0:2)] == 1, na.rm=TRUE) >= 1,
	bp_med_atc = apply(.[, paste0("20003-0.", 0:47)],1,function(x) any(x%in%med_atc_df$Coding[is_bp_med])),
	chol_med_atc = apply(.[, paste0("20003-0.", 0:47)],1,function(x) any(x%in%med_atc_df$Coding[is_chol_med])),
	bp_med = bp_med_sr | bp_med_atc, 
	chol_med = chol_med_sr | chol_med_atc
  ) 

bp_df <- basic_phenos_df %>%
  select(
    id = eid,
    dbp1 = `4079-0.0`, dbp2 = `4079-0.1`,
    sbp1 = `4080-0.0`, sbp2 = `4080-0.1`
  ) %>%
  inner_join(meds_df, by="id") %>%
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
  chol = 30690, glu = 30740, glu_date = 30741, 
  hba1c = 30750, hdl = 30760, ldl = 30780,
  tg = 30870, tg_date = 30871
)

biomarker_df <- left_join(blood_assays_df, basic_covar_df, by=c("eid" = "id")) %>% # merge in fasting_time
  select(
    id = eid,
    all_of(setNames(paste0(ukb_biomarker_fields, "-0.0"),
                   names(ukb_biomarker_fields))),
    fasting_time
  )

lipids_df <- biomarker_df %>%
  select(id, hdl, tg, ldl, glu_date, tg_date, fasting_time) %>%
  # Need fasting adjustment here (Note: glu.mg/dL = glu*18, chol.mg/dl = chol*38.67) for future reference
  inner_join(meds_df, by="id") %>%
  mutate(hdl = hdl * 38.67, # 1 mmol/L = 38.67 mg/dL for HDL
	 ldl = ldl * 38.67, # 1 mmol/L = 38.67 mg/dL for LDL
	 tg = tg * 88.57, # 1 mmol/L = 88.57 mg/dL for TG
	 fasting_status = ifelse(!is.na(fasting_time) & fasting_time >= 8, "fasted",
				 ifelse(!is.na(glu_date) & !is.na(tg_date) & glu_date == tg_date, 
					"pseudo-fasted", "non-fasted"))) %>%
  mutate(ldl_orig = ldl, # keeping original variables if needed later for descriptive statistics
	 hdl_orig = hdl, 
	 tg_orig = tg,
	 ldl = ifelse(chol_med, ldl / 0.7, ldl)) %>%
  mutate(across(c(hdl, tg), log)) %>%  # Confirm that log-transform should happen AFTER meds adjustment
  mutate(across(c(hdl, tg, ldl), winsorize)) %>%
  select(id, hdl, tg, ldl, hdl_orig, tg_orig, ldl_orig, chol_med, fasting_status)
  
outcome_df <- full_join(
  bp_df, 
  lipids_df,
  by="id"
)

# Final processing -------------------------------------------------------------

final_pheno_df <- covar_df %>%
  left_join(exposure_df, by="id") %>%
  left_join(outcome_df, by="id")

## winsorize CPD and PY by population and sex now that exposure and covar data are merged
final_pheno_df <- final_pheno_df %>%
  group_by(pop, sex) %>%
  mutate(across(c(cpd, py), winsorize))

pheno_filename <- paste0("ukb_phenos_", tag, ".csv")
write_csv(final_pheno_df, pheno_filename)
