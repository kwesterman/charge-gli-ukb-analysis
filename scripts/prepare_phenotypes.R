# Parse command-line argument(s) -----------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
tag <- ifelse(is.na(args[1]), "test", args[1]) # Usually should be the date in YYYYMMDD format

# Input filenames --------------------------------------------------------------

basic_phenos_file <- "../data/raw/ukb47880.csv"
blood_assays_file <- "../data/raw/ukb51002.csv"
urine_assays_file <- "../data/raw/ukb50582.csv"
pack_years_file <- "../data/raw/ukb668908.csv"
pa_met_file <- "../data/raw/ukb669852.csv"

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
  inner_join(pa_bridge_df, by="s")%>%
  filter(related == FALSE)

sample_exclusions <- readLines(sample_exclusion_file)  # Participants who revoked consent

valid_ids <- setdiff(pan_ancestry_df$id, sample_exclusions)  # Participants with genetics but not having revoked consent

pack_years_df <- fread(pack_years_file, data.table=FALSE) # Load pack years first to merge into basic pheno data
pa_met_df <- fread(pa_met_file, data.table=FALSE) # Load met first to merge into basic pheno data
basic_phenos_df <- fread(basic_phenos_file, data.table=FALSE) %>%
  left_join(pack_years_df, by="eid") %>%
  left_join(pa_met_df, by="eid") %>%
  rename('90016-0.0' = '90016-0.0.y') %>%
  filter(eid %in% valid_ids)
blood_assays_df <- fread(blood_assays_file, data.table=FALSE) %>%
  filter(eid %in% valid_ids)
urine_assays_df <- fread(urine_assays_file, data.table=FALSE) %>%
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
         light_vs_never = ifelse((!never & !light) | very_heavy, NA, ifelse(light, 1, 0)),
         heavy_vs_never = ifelse((!never & !heavy) | very_heavy, NA, ifelse(heavy, 1, 0)),
         curdrink = ifelse(very_heavy, NA, curdrink)) %>%
  select(id, curdrink, heavy_vs_light, light_vs_never, heavy_vs_never)

sleep_df <- basic_phenos_df %>%
  select(
    id = eid,
    age = `21022-0.0`, 
    sex = `31-0.0`,  # Coding shouldn't matter for regression adjustment
    sleep_dur = `1160-0.0`
  ) %>%
  mutate(sleep_dur = ifelse(sleep_dur %in% c(-1, -3), NA, sleep_dur)) %>%
  filter(!is.na(sleep_dur)) %>%
  inner_join(select(pan_ancestry_df, pop, id), by=c("id")) %>%
  filter(id %in% valid_ids) %>%
  filter(!duplicated(id))

sleep_df <- sleep_df %>%
  group_by(pop) %>%
  mutate(sleep_dur_resid = resid(lm(sleep_dur ~ sex * age))) %>% 
  mutate(stst = as.integer(sleep_dur_resid <= quantile(sleep_dur_resid, 0.2)),
         ltst = as.integer(sleep_dur_resid >= quantile(sleep_dur_resid, 0.8))) %>%
  ungroup() %>%
  select(id, sleep_dur_resid, stst, ltst)

source("depression_phenotyping.R")
depression_df <- readRDS("depression_intermediate_df.rds") %>% 
  inner_join(select(basic_phenos_df,`21022-0.0`, eid, `31-0.0`), by=c("f.eid"="eid")) %>% 
  inner_join(select(pan_ancestry_df, pop, id), by=c("f.eid"="id"))%>%
  filter(f.eid %in% valid_ids)%>%
  rename(age = '21022-0.0', 
         sex = '31-0.0') %>%
  mutate(
    dDEPR_MHQ = case_when(
      (has_mhq == 1) & (has_psychosis0 == 0) & 
        (PHQ_score >= 10) ~ 1,
      (has_mhq == 1) & (has_psychosis0 == 0) ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    dDEPR_nonMHQ = case_when(
      (has_mhq == 0) & (has_psychosis0 == 0) &
        (help_seeking_0 | sr_depression_0 | sr_antidepressant_usage_0 | 
           depression_hospital_icd10_0 | derived_single_md_0 | 
           derived_repeat_moderate_md_0 | derived_repeat_severe_md_0) ~ 1,
      (has_mhq == 0) & (has_psychosis0 == 0) ~ 0,
      TRUE ~ as.numeric(NA)
    ) 	
  ) %>%
  mutate(dDEPR = case_when(
    !is.na(dDEPR_MHQ) ~ dDEPR_MHQ,
    !is.na(dDEPR_nonMHQ) ~ dDEPR_nonMHQ,
    TRUE ~ as.numeric(NA)
  )) %>%  
  group_by(pop) %>%
  mutate(PHQ_score = winsorize(PHQ_score)) %>%
  nest() %>%
  rowwise() %>%
  mutate(
    data2 = list(calculate_qDEPR(data, "PHQ_score"))
  ) %>%
 unnest(data2) %>%
 ungroup() %>%
 select(id=f.eid, dDEPR_MHQ, dDEPR_nonMHQ, dDEPR, qDEPR)
 
  
pa_df <- basic_phenos_df %>%
  select(
    id = eid,
	sex = `31-0.0`,
    ipaq_met = `22040-0.0`, 
    walking_pleasure_frq = `971-0.0`, walking_pleasure_dur = `981-0.0`,
    strenuous_sport_frq = `991-0.0`, strenuous_sport_dur = `1001-0.0`, 
    heavy_diy_frq = `2624-0.0`, heavy_diy_dur = `2634-0.0`, 
    other_excercise_frq = `3637-0.0`, other_excercise_dur = `3647-0.0`, 
    contains("6164-0."),
    accel_avg = `90012-0.0`,
    accel_qc_good_wear_time = `90015-0.0`,
    accel_qc_good_calib = `90016-0.0`
  ) %>% 
  mutate(ipaq_met_p30_01 = ifelse((sex == 0 & ipaq_met > quantile(ipaq_met[sex==0],0.30,na.rm=T)) | 
								   (sex == 1 & ipaq_met > quantile(ipaq_met[sex==1],0.30,na.rm=T)), 1, 0)) %>%
  mutate(walking_pleasure_frq = ifelse(walking_pleasure_frq<0 | is.na(walking_pleasure_frq), 0 ,walking_pleasure_frq), # set neagatives and missings to 0
         walking_pleasure_dur = ifelse(walking_pleasure_dur<0 | is.na(walking_pleasure_dur), 0 ,walking_pleasure_dur),
         strenuous_sport_frq = ifelse(strenuous_sport_frq<0 | is.na(strenuous_sport_frq), 0 ,strenuous_sport_frq),
         strenuous_sport_dur = ifelse(strenuous_sport_dur<0 | is.na(strenuous_sport_dur), 0 ,strenuous_sport_dur),
         heavy_diy_frq = ifelse(heavy_diy_frq<0 | is.na(heavy_diy_frq), 0 ,heavy_diy_frq),
         heavy_diy_dur = ifelse(heavy_diy_dur<0 | is.na(heavy_diy_dur), 0 ,heavy_diy_dur),
         other_excercise_frq = ifelse(other_excercise_frq<0 | is.na(other_excercise_frq), 0 ,other_excercise_frq),
         other_excercise_dur = ifelse(other_excercise_dur<0 | is.na(other_excercise_dur), 0 ,other_excercise_dur),
  ) %>%
  mutate(rpaq_met = 3.3*walking_pleasure_frq*walking_pleasure_dur + 
           8*strenuous_sport_frq*strenuous_sport_dur + 
           4.5*heavy_diy_frq*heavy_diy_dur + 
           4.5*other_excercise_frq*other_excercise_dur
  ) %>%
  mutate(pa_any = apply(.[, paste0("6164-0.", 0:4)],1,function(x) any(x%in%c(1,2,3,5)))) %>% # exclude any 0's from above that indicated any PA in variable 6164
  mutate(rpaq_met = ifelse(pa_any & rpaq_met == 0, NA, rpaq_met)) %>% 
  mutate(rpaq_met_p30_01 = ifelse((sex == 0 & rpaq_met > quantile(rpaq_met[sex==0],0.30,na.rm=T)) | 
								   (sex == 1 & rpaq_met > quantile(rpaq_met[sex==1],0.30,na.rm=T)), 1, 0)) %>%
  mutate(accel_avg = ifelse(accel_qc_good_wear_time==0 | accel_qc_good_calib==0, NA, accel_avg)) %>%
  mutate(accel_avg_p30_01 = ifelse((sex == 0 & accel_avg > quantile(accel_avg[sex==0],0.30,na.rm=T)) | 
								   (sex == 1 & accel_avg > quantile(accel_avg[sex==1],0.30,na.rm=T)), 1, 0)) %>%
  select(id, rpaq_met, rpaq_met_p30_01)
  # select(id, ipaq_met, ipaq_met_p25_01, rpaq_met, rpaq_met_p25_01, accel_avg, accel_avg_p25_01)

source("glycemic_traits.R")
glycemic_df <-
  read_csv("ukb_glycemictraits.csv") %>%
  select(id, hba1c, glucose, fasting_status) %>% rename(glycemic_fasting_status=fasting_status)

source("diabetes.R") # this script doesn't run unless done interactively section by section, sorry
diabetes_df <-
  read_csv(paste0("ukb_diabetes.csv")) %>%
  select(id, diabetes_control, t2d_case, t1d_case)
			
education_df <- basic_phenos_df %>%
  select(id = eid)

exposure_df <- smoking_df %>%
  left_join(alcohol_df, by="id") %>%
  left_join(sleep_df, by="id") %>%
  left_join(depression_df, by="id") %>%
  left_join(pa_df, by="id") %>%
  left_join(education_df, by="id")%>%
  left_join(glycemic_df, by = "id") %>%
  left_join(diabetes_df, by = "id")

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

obesity_df <- basic_phenos_df %>%
  select(
    id = eid,
    bmi_wins = `21001-0.0`, waist = `48-0.0`,
    hip = `49-0.0`, pregnant = `3140-0.0`
  ) %>%
  mutate(
    whr_wins = waist / hip
  ) %>%
  mutate(across(c(bmi_wins, whr_wins), winsorize)) %>%
  filter(pregnant==0 | is.na(pregnant)) %>%
  select(id, bmi_wins, whr_wins)

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

outcome_df <- bp_df %>% 
  full_join(lipids_df,by="id") %>%
  full_join(obesity_df,by="id")

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
