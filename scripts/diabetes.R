library(tidyverse)
library(data.table)

filename <- "ukb_diabetes.csv"


# File names --------------------------------------------------------------


glycemic_file <- "../data/raw/ukb51002.csv"
basic_phenos_file <- "../data/raw/ukb47880.csv"
icd_file <- "../data/raw/ukb678359.csv"
add_icd_file <- "../data/raw/ukb671261.csv"

pan_ancestry_file <-
  "../data/raw/Files for retman/all_pops_non_eur_pruned_within_pop_pc_covs.tsv"
pa_bridge_file <- "../data/raw/ukb8343bridge31063.txt"
sample_exclusion_file <- "../data/raw/w8343_20220222.csv"


# Read in files -----------------------------------------------------------

pa_bridge_df <-
  fread(pa_bridge_file,
        col.names = c("id", "s"),
        data.table = FALSE)
pan_ancestry_df <- fread(pan_ancestry_file, data.table = FALSE) %>%
  inner_join(pa_bridge_df, by = "s") %>%
  filter(related == FALSE)

sample_exclusions <-
  readLines(sample_exclusion_file)  # Participants who revoked consent
valid_ids <-
  setdiff(pan_ancestry_df$id, sample_exclusions)  # Participants with genetics but not having revoked consent

diabetes_glycemic_df <- fread(glycemic_file, data.table = FALSE) %>%
  filter(eid %in% valid_ids)
##diabetes_basic_phenos_df <- fread(basic_phenos_file, data.table = FALSE) %>%
##  filter(eid %in% valid_ids)
icd_df <- fread(icd_file, data.table = FALSE) %>%
  rename(eid = f.eid) %>%
  filter(eid %in% valid_ids)
add_icd_df <- fread(add_icd_file, data.table = FALSE) %>%
  rename(eid = f.eid) %>%
  filter(eid %in% valid_ids)

# Hba1c -------------------------------------------------------------------


a1c <- diabetes_glycemic_df %>%
  select(id = eid,
         hba1c0 = `30750-0.0`,
         hba1c1 = `30750-1.0`) %>%
  mutate(
    hba1c.30750.IFCC.max = pmax(hba1c0, hba1c1, na.rm = TRUE),
    hba1c.30750.NGSP.max = hba1c.30750.IFCC.max / 10.929 + 2.15
  )



# Diabetes ----------------------------------------------------------------

feids <-
  c("eid", "20002","20003","4041","20009","2976","6177","6153","31","34","52",
    "21022","21003","53","54","3140","50","49","48","23098","23099","23100",
    "23101","23102","23105","21002","21001","4080","93","4079","94","30515",
    "30510","30505","30500","2986","2976","2443","5901")

diabetes_df <- basic_phenos_df %>%
  select(contains(feids))

feids <- c("eid", "20107", "20110", "20111", "5901")

diabetes_fam_hist <- icd_df %>%
  select(contains(feids))

colnames(diabetes_fam_hist) <-
  gsub("f\\.", "", colnames(diabetes_fam_hist))
colnames(diabetes_fam_hist) <-
  sub("\\.", "-", colnames(diabetes_fam_hist))

diabetes_df <- diabetes_df %>%
  left_join(diabetes_fam_hist, by = "eid")

diabetes_df$father_dm <- as.integer(Reduce('|', lapply(diabetes_df[paste('20107-', c(0:3),'.', c(0:9), sep = "")], '%in%', c(9)))) 
table(diabetes_df$father_dm, useNA = "always")

diabetes_df$mother_dm <- as.integer(Reduce('|', lapply(diabetes_df[paste('20110-', c(0:3),'.', c(0:9), sep = "")], '%in%', c(9)))) 
table(diabetes_df$mother_dm, useNA = "always")

diabetes_df$sibling_dm <- as.integer(Reduce('|', lapply(diabetes_df[paste('20111-', c(0:3),'.', c(0:9), sep = "")], '%in%', c(9)))) 
table(diabetes_df$sibling_dm, useNA = "always")


# ICD codes --------------------------------------------------------------------
feids <- c("eid", "41203", "41205")

icd_codes <- icd_df %>%
  select(contains(feids)) %>% rename(id = eid)

feids <- c("eid", "eid", "41202", "41204")
icd_10 <- basic_phenos_df %>%
  select(contains(feids)) %>%
  left_join(add_icd_df %>%  select(contains(feids)),
            by = c("eid" = "eid"))

icd_codes <- icd_codes %>% left_join(icd_10, by = c("id" = "eid"))

colnames(icd_codes) <- gsub("f\\.", "", colnames(icd_codes))
colnames(icd_codes) <- sub("\\.0\\.", "-0\\.", colnames(icd_codes))

codes <- c("E100","E101","E102","E103","E104","E105","E107","E108","E109")
icd_codes$E10_main <- as.integer(Reduce('|', lapply(icd_codes[paste('41202-0.', c(0:74), sep = "")], '%in%', codes))) 
table(icd_codes$E10_main, useNA = "always")
icd_codes$E10_secondary <- as.integer(Reduce('|', lapply(icd_codes[paste('41204-0.', c(0:187), sep = "")], '%in%', codes))) 
table(icd_codes$E10_secondary, useNA = "always")

codes <-  c("2500","2501","2502","2503","2504","2505","2507","2508","2509")
icd_codes$icd9_main <- as.integer(Reduce('|', lapply(icd_codes[paste('41203-0.', c(0:27), sep = "")], '%in%', codes)))
table(icd_codes$icd9_main, useNA = "always")
icd_codes$icd9_secondary <- as.integer(Reduce('|', lapply(icd_codes[paste('41205-0.', c(0:29), sep = "")], '%in%', codes))) 
table(icd_codes$icd9_secondary, useNA = "always")

codes <- c("E080","E081","E082","E083","E084","E085","E087","E088","E089","E090","E091","E092",
           "E093","E094","E095","E097","E098","E099","E100","E101","E102","E103","E104","E105","E107",
           "E108","E109","E110","E111","E112","E113","E114","E115","E117","E118","E119","E130","E131",
           "E132","E133","E134","E135","E137","E138","E139")
icd_codes$icd10_main <- as.integer(Reduce('|', lapply(icd_codes[paste('41202-0.', c(0:65), sep = "")], '%in%', codes))) 
table(icd_codes$icd10_main, useNA = "always")
icd_codes$icd10_secondary <- as.integer(Reduce('|', lapply(icd_codes[paste('41204-0.', c(0:65), sep = "")], '%in%', codes))) 
table(icd_codes$icd10_secondary, useNA = "always")

df <- diabetes_df

# DM ----------------------------------------------------------------------

df$dm_any_sr_ni.0 <- as.integer(Reduce('|', lapply(df[paste('20002-0.', c(0:33), sep = "")], '%in%', c(1220)))) 
table(df$dm_any_sr_ni.0, useNA = "always")

df$dm_any_sr_ni.1 <- as.integer(Reduce('|', lapply(df[paste('20002-1.', c(0:33), sep = "")], '%in%', c(1220)))) 
table(df$dm_any_sr_ni.1, useNA = "always")

df$dm_any_sr_ni.2 <- as.integer(Reduce('|', lapply(df[paste('20002-2.', c(0:33), sep = "")], '%in%', c(1220)))) 
table(df$dm_any_sr_ni.2, useNA="always")

df$dm_any_sr_ni.3 <- as.integer(Reduce('|', lapply(df[paste('20002-3.', c(0:33), sep = "")], '%in%', c(1220)))) 
table(df$dm_any_sr_ni.3, useNA="always")

df$dm_any_sr_ni =ifelse(df$dm_any_sr_ni.0 ==1 |df$dm_any_sr_ni.1 ==1 | df$dm_any_sr_ni.2 ==1 |df$dm_any_sr_ni.3 ==1,
                        1, 0)
table(df$dm_any_sr_ni, useNA="always")

# GDM ---------------------------------------------------------------------

df$dm_gdm_sr_ni.0 <- as.integer(Reduce('|', lapply(df[paste('20002-0.', c(0:33), sep = "")], '%in%', c(1221)))) 
table(df$dm_gdm_sr_ni.0, useNA="always")
df$dm_gdm_sr_ni.1 <- as.integer(Reduce('|', lapply(df[paste('20002-1.', c(0:33), sep = "")], '%in%', c(1221)))) 
table(df$dm_gdm_sr_ni.1, useNA="always")
df$dm_gdm_sr_ni.2 <- as.integer(Reduce('|', lapply(df[paste('20002-2.', c(0:33), sep = "")], '%in%', c(1221))))  
table(df$dm_gdm_sr_ni.2, useNA="always")
df$dm_gdm_sr_ni.3 <- as.integer(Reduce('|', lapply(df[paste('20002-3.', c(0:33), sep = "")], '%in%', c(1221)))) 
table(df$dm_gdm_sr_ni.3, useNA="always")

df$dm_gdm_sr_ni = ifelse(df$dm_gdm_sr_ni.0 ==1 |df$dm_gdm_sr_ni.1 ==1 | df$dm_gdm_sr_ni.2 ==1 | df$dm_gdm_sr_ni.3 ==1,
                         1, 0)
table(df$dm_gdm_sr_ni, useNA="always")


# T1D ---------------------------------------------------------------------

df$dm_t1dm_sr_ni.0 <- as.integer(Reduce('|', lapply(df[paste('20002-0.', c(0:33), sep = "")], '%in%', c(1222)))) 
table(df$dm_t1dm_sr_ni.0, useNA="always")
df$dm_t1dm_sr_ni.1 <- as.integer(Reduce('|', lapply(df[paste('20002-1.', c(0:33), sep = "")], '%in%', c(1222)))) 
table(df$dm_t1dm_sr_ni.1, useNA="always")
df$dm_t1dm_sr_ni.2 <- as.integer(Reduce('|', lapply(df[paste('20002-2.', c(0:33), sep = "")], '%in%', c(1222)))) 
table(df$dm_t1dm_sr_ni.2, useNA="always")
df$dm_t1dm_sr_ni.3 <- as.integer(Reduce('|', lapply(df[paste('20002-3.', c(0:33), sep = "")], '%in%', c(1222)))) 
table(df$dm_t1dm_sr_ni.3, useNA="always")

df$dm_t1dm_sr_ni = ifelse(df$dm_t1dm_sr_ni.0 ==1 |df$dm_t1dm_sr_ni.1 ==1 |df$dm_t1dm_sr_ni.2 ==1 |df$dm_t1dm_sr_ni.3 ==1 ,
                          1, 0)
table(df$dm_t1dm_sr_ni, useNA="always")



# T2D ---------------------------------------------------------------------

df$dm_t2dm_sr_ni.0 <- as.integer(Reduce('|', lapply(df[paste('20002-0.', c(0:33), sep = "")], '%in%', c(1223)))) 
table(df$dm_t2dm_sr_ni.0, useNA="always")
df$dm_t2dm_sr_ni.1 <- as.integer(Reduce('|', lapply(df[paste('20002-1.', c(0:33), sep = "")], '%in%', c(1223)))) 
table(df$dm_t2dm_sr_ni.1, useNA="always")
df$dm_t2dm_sr_ni.2 <- as.integer(Reduce('|', lapply(df[paste('20002-2.', c(0:33), sep = "")], '%in%', c(1223)))) 
table(df$dm_t2dm_sr_ni.2, useNA="always")
df$dm_t2dm_sr_ni.3 <- as.integer(Reduce('|', lapply(df[paste('20002-3.', c(0:33), sep = "")], '%in%', c(1223)))) 
table(df$dm_t2dm_sr_ni.3, useNA="always")

df$dm_t2dm_sr_ni = ifelse(df$dm_t2dm_sr_ni.0 ==1 |df$dm_t2dm_sr_ni.1 ==1 |df$dm_t2dm_sr_ni.2 ==1 |df$dm_t2dm_sr_ni.3 ==1 ,
                          1, 0)
table(df$dm_t2dm_sr_ni, useNA="always")


# Non-specific DM ---------------------------------------------------------

## Replacing 10 people with both T1 and T2 diagnoses to non-specific
## these people now have t1t2dm + dm, but have 0's for t1dm and t2dm specifically
df$t1t2dm_sr_ni <- ifelse((df$dm_t1dm_sr_ni == 1 & df$dm_t2dm_sr_ni == 1), 1, 0)
table(df$t1t2dm_sr_ni, useNA="always")

table(df$dm_any_sr_ni, useNA="always")
df$dm_any_sr_ni <- ifelse(df$t1t2dm_sr_ni==1, 1, df$dm_any_sr_ni)
table(df$dm_any_sr_ni, useNA="always")

table(df$dm_t1dm_sr_ni, useNA="always")
df$dm_t1dm_sr_ni <- ifelse(df$t1t2dm_sr_ni==1, 0, df$dm_t1dm_sr_ni)
table(df$dm_t1dm_sr_ni, useNA="always")

table(df$dm_t2dm_sr_ni, useNA="always")
df$dm_t2dm_sr_ni <- ifelse(df$t1t2dm_sr_ni==1, 0, df$dm_t2dm_sr_ni)
table(df$dm_t2dm_sr_ni, useNA="always")


# TS GDM ------------------------------------------------------------------

table(df$`4041-0.0`, useNA="always")

df$dm_gdm_sr_ts.0 <- ifelse(df$`4041-0.0` != 1 | is.na(df$`4041-0.0`), 0,1)
table(df$dm_gdm_sr_ts.0, useNA="always")
df$dm_gdm_sr_ts.1 <- ifelse(df$`4041-1.0` != 1 | is.na(df$`4041-1.0`), 0,1)
table(df$dm_gdm_sr_ts.1, useNA="always")
df$dm_gdm_sr_ts.2 <- ifelse(df$`4041-2.0` != 1 | is.na(df$`4041-2.0`), 0,1)
table(df$dm_gdm_sr_ts.2, useNA="always")
df$dm_gdm_sr_ts.3 <- ifelse(df$`4041-3.0` != 1| is.na(df$`4041-3.0`), 0,1)
table(df$dm_gdm_sr_ts.3, useNA="always")

df$dm_gdm_sr_ts <- ifelse(df$dm_gdm_sr_ts.0 == 1 | df$dm_gdm_sr_ts.0 == 1 |df$dm_gdm_sr_ts.0 == 1, 
                          1,0)
table(df$dm_gdm_sr_ts, useNA="always")	

table(df$dm_gdm_sr_ts, df$dm_gdm_sr_ni)

# Age at diabetes NI ---------------------------------------------------------
df$age_dm_any_dx_sr.0 <- -1
for (i in 0:33) {df$age_dm_any_dx_sr.0 <- ifelse(is.na(df[[paste('20002-0.', i, sep="")]]) == T | 
                                                   df[[paste('20002-0.', i, sep="")]] != 1220, 
                                                 df$age_dm_any_dx_sr.0, df[[paste('20009-0.', i, sep="")]])}
df$age_dm_any_dx_sr.0 <- ifelse(df$age_dm_any_dx_sr.0<0 | is.na(df$age_dm_any_dx_sr.0), NA, df$age_dm_any_dx_sr.0)
table(df$age_dm_any_dx_sr.0, useNA="always")

df$age_dm_gdm_dx_sr.0 <- -1
for (i in 0:33) {df$age_dm_gdm_dx_sr.0 <- ifelse(is.na(df[[paste('20002-0.', i, sep="")]])| 
                                                   df[[paste('20002-0.', i, sep="")]] != 1221, 
                                                 df$age_dm_gdm_dx_sr.0, df[[paste('20009-0.', i, sep="")]])}
df$age_dm_gdm_dx_sr.0 <- ifelse(df$age_dm_gdm_dx_sr.0<0 | is.na(df$age_dm_gdm_dx_sr.0), NA, df$age_dm_gdm_dx_sr.0)
table(df$age_dm_gdm_dx_sr.0, useNA="always")

df$age_dm_t1dm_dx_sr.0 <- -1
for (i in 0:33) {df$age_dm_t1dm_dx_sr.0 <- ifelse(is.na(df[[paste('20002-0.', i, sep="")]])| 
                                                    df[[paste('20002-0.', i, sep="")]] != 1222, 
                                                  df$age_dm_t1dm_dx_sr.0, df[[paste('20009-0.', i, sep="")]])}
df$age_dm_t1dm_dx_sr.0 <- ifelse(df$age_dm_t1dm_dx_sr.0 <0 | is.na(df$age_dm_t1dm_dx_sr.0) , NA, df$age_dm_t1dm_dx_sr.0)
table(df$age_dm_t1dm_dx_sr.0, useNA="always")

df$age_dm_t2dm_dx_sr.0 <- -1
for (i in 0:33) {df$age_dm_t2dm_dx_sr.0 <- ifelse(is.na(df[[paste('20002-0.', i, sep="")]])|
                                                    df[[paste('20002-0.', i, sep="")]] != 1223, 
                                                  df$age_dm_t2dm_dx_sr.0, df[[paste('20009-0.', i, sep="")]])}
df$age_dm_t2dm_dx_sr.0 <- ifelse(df$age_dm_t2dm_dx_sr.0 <0 | is.na(df$age_dm_t2dm_dx_sr.0), NA, df$age_dm_t2dm_dx_sr.0)
table(df$age_dm_t2dm_dx_sr.0, useNA="always")

df$age_dm_any_dx_sr.1 <- -1
for (i in 0:33) {df$age_dm_any_dx_sr.1 <- ifelse(is.na(df[[paste('20002-1.', i, sep="")]]) |
                                                   df[[paste('20002-1.', i, sep="")]] != 1220, 
                                                 df$age_dm_any_dx_sr.1, df[[paste('20009-1.', i, sep="")]])}
df$age_dm_any_dx_sr.1 <- ifelse(df$age_dm_any_dx_sr.1<0 | is.na(df$age_dm_any_dx_sr.1), NA, df$age_dm_any_dx_sr.1)
table(df$age_dm_any_dx_sr.1, useNA="always")

df$age_dm_gdm_dx_sr.1 <- -1
for (i in 0:33) {df$age_dm_gdm_dx_sr.1 <- ifelse(is.na(df[[paste('20002-1.', i, sep="")]])|
                                                   df[[paste('20002-1.', i, sep="")]] != 1221, 
                                                 df$age_dm_gdm_dx_sr.1, df[[paste('20009-1.', i, sep="")]])}
df$age_dm_gdm_dx_sr.1 <- ifelse(df$age_dm_gdm_dx_sr.1<0 | is.na(df$age_dm_gdm_dx_sr.1) , NA, df$age_dm_gdm_dx_sr.1)
table(df$age_dm_gdm_dx_sr.1, useNA="always")
df$age_dm_t1dm_dx_sr.1 <- -1
for (i in 0:33) {df$age_dm_t1dm_dx_sr.1 <- ifelse(is.na(df[[paste('20002-1.', i, sep="")]])|
                                                    df[[paste('20002-1.', i, sep="")]] != 1222, 
                                                  df$age_dm_t1dm_dx_sr.1, df[[paste('20009-1.', i, sep="")]])}
df$age_dm_t1dm_dx_sr.1 <- ifelse(df$age_dm_t1dm_dx_sr.1 <0 | is.na(df$age_dm_t1dm_dx_sr.1), NA, df$age_dm_t1dm_dx_sr.1)
table(df$age_dm_t1dm_dx_sr.1, useNA="always")
df$age_dm_t2dm_dx_sr.1 <- -1
for (i in 0:33) {df$age_dm_t2dm_dx_sr.1 <- ifelse(is.na(df[[paste('20002-1.', i, sep="")]])|
                                                    df[[paste('20002-1.', i, sep="")]] != 1223, 
                                                  df$age_dm_t2dm_dx_sr.1, df[[paste('20009-1.', i, sep="")]])}
df$age_dm_t2dm_dx_sr.1 <- ifelse(df$age_dm_t2dm_dx_sr.1 <0 | is.na(df$age_dm_t2dm_dx_sr.1), NA, df$age_dm_t2dm_dx_sr.1)
table(df$age_dm_t2dm_dx_sr.1, useNA="always")
df$age_dm_any_dx_sr.2 <- -1
for(i in 0:33) {df$age_dm_any_dx_sr.2 <- ifelse(is.na(df[[paste('20002-2.', i, sep="")]]) | 
                                                   df[[paste('20002-2.', i, sep="")]] != 1220, 
                                                 df$age_dm_any_dx_sr.2, df[[paste('20009-2.', i, sep="")]])}

df$age_dm_any_dx_sr.2 <- ifelse(df$age_dm_any_dx_sr.2<0 | is.na(df$age_dm_any_dx_sr.2), NA, df$age_dm_any_dx_sr.2)
table(df$age_dm_any_dx_sr.2, useNA="always")

df$age_dm_gdm_dx_sr.2 <- -1
for(i in 0:33) {df$age_dm_gdm_dx_sr.2 <- ifelse(is.na(df[[paste('20002-2.', i, sep="")]])|
                                                  df[[paste('20002-2.', i, sep="")]] != 1221, 
                                                df$age_dm_any_dx_sr.2, df[[paste('20009-2.', i, sep="")]])}
df$age_dm_gdm_dx_sr.2 <- ifelse(df$age_dm_gdm_dx_sr.2<0|is.na(df$age_dm_gdm_dx_sr.2),NA,df$age_dm_gdm_dx_sr.2)
table(df$age_dm_gdm_dx_sr.2, useNA="always")

df$age_dm_t1dm_dx_sr.2 <- -1
for(i in 0:33){df$age_dm_t1dm_dx_sr.2<-ifelse(is.na(df[[paste('20002-2.', i, sep="")]])|df[[paste('20002-2.', i, sep="")]] != 1222, df$age_dm_t1dm_dx_sr.2, df[[paste('20009-2.', i, sep="")]])}
df$age_dm_t1dm_dx_sr.2 <- ifelse(df$age_dm_t1dm_dx_sr.2 <0 | is.na(df$age_dm_t1dm_dx_sr.2) , NA, df$age_dm_t1dm_dx_sr.2)
table(df$age_dm_t1dm_dx_sr.2, useNA="always")

df$age_dm_t2dm_dx_sr.2 <- -1
for(i in 0:33){df$age_dm_t2dm_dx_sr.2<-ifelse(is.na(df[[paste('20002-2.', i, sep="")]])|
                                                    df[[paste('20002-2.', i, sep="")]] != 1223,
                                                  df$age_dm_t2dm_dx_sr.2, df[[paste('20009-2.', i, sep="")]])}
df$age_dm_t2dm_dx_sr.2<-ifelse(df$age_dm_t2dm_dx_sr.2 <0 | is.na(df$age_dm_t2dm_dx_sr.2) , NA, df$age_dm_t2dm_dx_sr.2)
table(df$age_dm_t2dm_dx_sr.2, useNA="always")

df$age_dm_any_dx_sr.3<- -1
for (i in 0:33) {df$age_dm_any_dx_sr.3<-ifelse(is.na(df[[paste('20002-3.', i, sep="")]])|df[[paste('20002-3.', i, sep="")]] != 1220, df$age_dm_any_dx_sr.3, df[[paste('20009-3.', i, sep="")]])}
df$age_dm_any_dx_sr.3 <- ifelse(df$age_dm_any_dx_sr.3<0 | is.na(df$age_dm_any_dx_sr.3) , NA, df$age_dm_any_dx_sr.3)
table(df$age_dm_any_dx_sr.3, useNA="always")
df$age_dm_gdm_dx_sr.3 <- -1
for (i in 0:33) {df$age_dm_gdm_dx_sr.3 <- ifelse(is.na(df[[paste('20002-3.', i, sep="")]])|df[[paste('20002-3.', i, sep="")]] != 1221, df$age_dm_gdm_dx_sr.3, df[[paste('20009-3.', i, sep="")]])}
df$age_dm_gdm_dx_sr.3 <- ifelse(df$age_dm_gdm_dx_sr.3<0 | is.na(df$age_dm_gdm_dx_sr.3) , NA, df$age_dm_gdm_dx_sr.3)
table(df$age_dm_gdm_dx_sr.3, useNA="always")
df$age_dm_t1dm_dx_sr.3 <- -1
for (i in 0:33) {df$age_dm_t1dm_dx_sr.3 <- ifelse(is.na(df[[paste('20002-3.', i, sep="")]])|df[[paste('20002-3.', i, sep="")]] != 1222, df$age_dm_t1dm_dx_sr.3, df[[paste('20009-3.', i, sep="")]])}
df$age_dm_t1dm_dx_sr.3 <- ifelse(df$age_dm_t1dm_dx_sr.3 <0 | is.na(df$age_dm_t1dm_dx_sr.3), NA, df$age_dm_t1dm_dx_sr.3)
table(df$age_dm_t1dm_dx_sr.3, useNA="always")
df$age_dm_t2dm_dx_sr.3 <- -1
for (i in 0:33) {df$age_dm_t2dm_dx_sr.3 <- ifelse(is.na(df[[paste('20002-3.', i, sep="")]])|df[[paste('20002-3.', i, sep="")]] != 1223, df$age_dm_t2dm_dx_sr.3, df[[paste('20009-3.', i, sep="")]])}
df$age_dm_t2dm_dx_sr.3 <- ifelse(df$age_dm_t2dm_dx_sr.3 <0 | is.na(df$age_dm_t2dm_dx_sr.3) , NA, df$age_dm_t2dm_dx_sr.3)
table(df$age_dm_t2dm_dx_sr.3, useNA="always")


### average diagnosis age - SR, NI. any, t1dm, t2dm, gdm -- 20009 associated with 20002
df$age_dm_any_dx_sr = rowMeans(df[,c("age_dm_any_dx_sr.0", "age_dm_any_dx_sr.1", "age_dm_any_dx_sr.2", "age_dm_any_dx_sr.3")],na.rm = TRUE)
table(df$age_dm_any_dx_sr, useNA="always")
df$age_dm_any_dx_sr[is.nan(df$age_dm_any_dx_sr)] <-NA
table(df$age_dm_any_dx_sr, useNA="always")

df$age_dm_gdm_dx_sr = rowMeans(df[,c("age_dm_gdm_dx_sr.0", "age_dm_gdm_dx_sr.1", "age_dm_gdm_dx_sr.2", "age_dm_gdm_dx_sr.3")],na.rm = TRUE)
table(df$age_dm_gdm_dx_sr, useNA="always")
df$age_dm_gdm_dx_sr[is.nan(df$age_dm_gdm_dx_sr)] <-NA
table(df$age_dm_gdm_dx_sr, useNA="always")

df$age_dm_t1dm_dx_sr = rowMeans(df[,c("age_dm_t1dm_dx_sr.0", "age_dm_t1dm_dx_sr.1", "age_dm_t1dm_dx_sr.2", "age_dm_t1dm_dx_sr.3")],na.rm = TRUE)
table(df$age_dm_t1dm_dx_sr, useNA="always")
df$age_dm_t1dm_dx_sr[is.nan(df$age_dm_t1dm_dx_sr)] <-NA
table(df$age_dm_t1dm_dx_sr, useNA="always")

df$age_dm_t2dm_dx_sr = rowMeans(df[,c("age_dm_t2dm_dx_sr.0", "age_dm_t2dm_dx_sr.1", "age_dm_t2dm_dx_sr.2", "age_dm_t2dm_dx_sr.3")],na.rm = TRUE)
table(df$age_dm_t2dm_dx_sr, useNA="always")
df$age_dm_t2dm_dx_sr [is.nan(df$age_dm_t2dm_dx_sr)] <- NA
table(df$age_dm_t2dm_dx_sr, useNA="always")


# Age variable combined TS and NI -----------------------------------------
#2976-0.0	Age diabetes diagnosed == HSD code calls this "T2D_DX_AGE"

##### set -1 and -3 to NA
df$`2976-0.0` = ifelse(df$`2976-0.0` < 0, NA, df$`2976-0.0`)
table(df$`2976-0.0`, useNA="always")
df$`2976-1.0` = ifelse(df$`2976-1.0` < 0, NA, df$`2976-1.0`)
table(df$`2976-1.0`, useNA="always")
df$`2976-2.0` = ifelse(df$`2976-2.0` < 0, NA, df$`2976-2.0`)
table(df$`2976-2.0`, useNA="always")
df$`2976-3.0` = ifelse(df$`2976-3.0` < 0, NA, df$`2976-3.0`)
table(df$`2976-3.0`, useNA="always")

#average Age diabetes diagnosed

df$`2976` = rowMeans(df[,c("2976-0.0", "2976-1.0", "2976-2.0", "2976-3.0")],na.rm = TRUE)
df$`2976`[is.nan(df$`2976`)] <- NA
table(df$`2976`, useNA="always")


### if you have T1DM -> age at T1DM
### if you have T2DM -> age at T2DM
### if age at T1DM, T2DM, and any DM is NA, AND you do not have gestational diabetes (as recorded on TS) -> TS age at dx (2976)
### this is because "Field 2976 was collected from men who indicated that a doctor had told them 
### they have diabetes, as defined by their answers to Field 2443 and all women except those who 
### indicated they had diabetes only during pregnancy, as defined by their answers to Field 2443
df$agedm_ts_or_ni <- ifelse(df$dm_t2dm_sr_ni == 1, df$age_dm_t2dm_dx_sr, 
                            ifelse(df$dm_t1dm_sr_ni == 1, df$age_dm_t1dm_dx_sr, 
                                   ifelse(is.na(df$age_dm_any_dx_sr)==T & is.na(df$age_dm_t1dm_dx_sr)==T  & is.na(df$age_dm_t2dm_dx_sr)==T & df$dm_gdm_sr_ts==0, df$`2976`, df$age_dm_any_dx_sr)))
table(df$agedm_ts_or_ni, useNA="always")
df$agedm_ts_or_ni <- ifelse(is.na(df$`2976`) ==T & is.na(df$age_dm_any_dx_sr)==T & is.na(df$age_dm_t1dm_dx_sr)==T & is.na(df$age_dm_t2dm_dx_sr)==T, NA, df$agedm_ts_or_ni)
table(df$agedm_ts_or_ni, useNA="always")	
df$agedm_ts_or_ni <- ifelse(df$agedm_ts_or_ni < 0, NA, df$agedm_ts_or_ni)
table(df$agedm_ts_or_ni, useNA="always")


# Insulin -----------------------------------------------------------------

### this is where field 6153 plays in - the field for insulin and other medications for women!
### so use field 6177 and 6153 for both men and women	
#Creating Self report diabetes medication: Touchscreen (TS):
df$dm_insulin_sr_ts <- 0 #use this single variable for self reported (via touchscreen) insulin use for both sexes (M and F)

df$dm_insulin_sr_ts.0 <- 0
for (i in 0:2){df$dm_insulin_sr_ts.0[df[[paste('6177-0.', i, sep="")]] == 3] <- 1 } 
table(df$dm_insulin_sr_ts.0, useNA="always")
for (i in 0:3){df$dm_insulin_sr_ts.0[df[[paste('6153-0.', i, sep="")]] == 3] <- 1 } 
table(df$dm_insulin_sr_ts.0, useNA="always")

df$dm_insulin_sr_ts.1 <- 0
for (i in 0:2){df$dm_insulin_sr_ts.1[df[[paste('6177-1.', i, sep="")]] == 3] <- 1 }          
table(df$dm_insulin_sr_ts.1, useNA="always")
for (i in 0:3){df$dm_insulin_sr_ts.1[df[[paste('6153-1.', i, sep="")]] == 3] <- 1 } 
table(df$dm_insulin_sr_ts.1, useNA="always")

df$dm_insulin_sr_ts.2 <- 0
for (i in 0:2){df$dm_insulin_sr_ts.2[df[[paste('6177-2.', i, sep="")]] == 3] <- 1 }          
table(df$dm_insulin_sr_ts.2, useNA="always")
for (i in 0:3){df$dm_insulin_sr_ts.2[df[[paste('6153-2.', i, sep="")]] == 3] <- 1 } 
table(df$dm_insulin_sr_ts.2, useNA="always")

df$dm_insulin_sr_ts.3 <- 0
for (i in 0:2){df$dm_insulin_sr_ts.3[df[[paste('6177-3.', i, sep="")]] == 3] <- 1 }          
table(df$dm_insulin_sr_ts.3, useNA="always")
for (i in 0:3){df$dm_insulin_sr_ts.3[df[[paste('6153-3.', i, sep="")]] == 3] <- 1 } 
table(df$dm_insulin_sr_ts.3, useNA="always")

df$dm_insulin_sr_ts = ifelse(df$dm_insulin_sr_ts.0 == 1 |
                               df$dm_insulin_sr_ts.1 == 1 |
                               df$dm_insulin_sr_ts.2 == 1 |
                               df$dm_insulin_sr_ts.3 == 1, 1, 0)
table(df$dm_insulin_sr_ts, useNA="always")



# Creating Medication Code NI ---------------------------------------------
# 3 categories of medications: insulin, metformin, non-metformin OAD
df$meds_insulin_sr_ni.0 <- 0
for (i in 0:47){ df$meds_insulin_sr_ni.0[df[[paste('20003-0.', i, sep="")]] == 1140883066] <- 1 }          
df$meds_metformin_sr_ni.0 <- 0
table(df$meds_insulin_sr_ni.0, useNA="always")

v <- c(1140884600, 1140874686, 1141189090)
for (i in 0:47){ df$meds_metformin_sr_ni.0[df[[paste('20003-0.', i, sep="")]]  %in% v] <- 1 }          
table(df$meds_metformin_sr_ni.0, useNA="always")

df$meds_nonmet_oad_sr_ni.0 <- 0
v <- c(1140874718, 1140874744, 1140874746, 1141152590, 1141156984, 1140874646, 1141157284, 1140874652, 1140874674, 1140874728, 1140868902, 1140868908, 1140857508, 1141173882, 1141173786, 1141168660, 1141171646, 1141171652, 1141153254, 1141177600, 1141177606)
for (i in 0:47){ df$meds_nonmet_oad_sr_ni.0[df[[paste('20003-0.', i, sep="")]] %in% v] <- 1 }
table(df$meds_nonmet_oad_sr_ni.0, useNA="always")

df$meds_insulin_sr_ni.1 <- 0
for (i in 0:47){ df$meds_insulin_sr_ni.1[df[[paste('20003-1.', i, sep="")]] == 1140883066] <- 1 }          
df$meds_metformin_sr_ni.1 <- 0
table(df$meds_insulin_sr_ni.1, useNA="always")

v <- c(1140884600, 1140874686, 1141189090)
for (i in 0:47){ df$meds_metformin_sr_ni.1[df[[paste('20003-1.', i, sep="")]]  %in% v] <- 1 }          
table(df$meds_metformin_sr_ni.1, useNA="always")

df$meds_nonmet_oad_sr_ni.1 <- 0
v <- c(1140874718, 1140874744, 1140874746, 1141152590, 1141156984, 1140874646, 1141157284, 1140874652, 1140874674, 1140874728, 1140868902, 1140868908, 1140857508, 1141173882, 1141173786, 1141168660, 1141171646, 1141171652, 1141153254, 1141177600, 1141177606)
for (i in 0:47){ df$meds_nonmet_oad_sr_ni.1[df[[paste('20003-1.', i, sep="")]] %in% v] <- 1 }
table(df$meds_nonmet_oad_sr_ni.1, useNA="always")

df$meds_insulin_sr_ni.2 <- 0
for (i in 0:47){ df$meds_insulin_sr_ni.2[df[[paste('20003-2.', i, sep="")]] == 1140883066] <- 1 }          

df$meds_metformin_sr_ni.2 <- 0
table(df$meds_insulin_sr_ni.2, useNA="always")

v <- c(1140884600, 1140874686, 1141189090)
for (i in 0:47){ df$meds_metformin_sr_ni.2[df[[paste('20003-2.', i, sep="")]]  %in% v] <- 1 }          
table(df$meds_metformin_sr_ni.2, useNA="always")

df$meds_nonmet_oad_sr_ni.2 <- 0
v <- c(1140874718, 1140874744, 1140874746, 1141152590, 1141156984, 1140874646, 1141157284, 1140874652, 1140874674, 1140874728, 1140868902, 1140868908, 1140857508, 1141173882, 1141173786, 1141168660, 1141171646, 1141171652, 1141153254, 1141177600, 1141177606)
for (i in 0:47){ df$meds_nonmet_oad_sr_ni.2[df[[paste('20003-2.', i, sep="")]] %in% v] <- 1 }
table(df$meds_nonmet_oad_sr_ni.2, useNA="always")


df$meds_insulin_sr_ni.3 <- 0
for (i in 0:47){ df$meds_insulin_sr_ni.3[df[[paste('20003-3.', i, sep="")]] == 1140883066] <- 1 }          
df$meds_metformin_sr_ni.3 <- 0
table(df$meds_insulin_sr_ni.3, useNA="always")

v <- c(1140884600, 1140874686, 1141189090)
for (i in 0:47){ df$meds_metformin_sr_ni.3[df[[paste('20003-3.', i, sep="")]]  %in% v] <- 1 }          
table(df$meds_metformin_sr_ni.3, useNA="always")

df$meds_nonmet_oad_sr_ni.3 <- 0
v <- c(1140874718, 1140874744, 1140874746, 1141152590, 1141156984, 1140874646, 1141157284, 1140874652, 1140874674, 1140874728, 1140868902, 1140868908, 1140857508, 1141173882, 1141173786, 1141168660, 1141171646, 1141171652, 1141153254, 1141177600, 1141177606)
for (i in 0:47){ df$meds_nonmet_oad_sr_ni.3[df[[paste('20003-3.', i, sep="")]] %in% v] <- 1 }
table(df$meds_nonmet_oad_sr_ni.3, useNA="always")


df$meds_insulin_sr_ni = ifelse(df$meds_insulin_sr_ni.0 == 1 |
                                 df$meds_insulin_sr_ni.1 == 1 |
                                 df$meds_insulin_sr_ni.2 == 1 |
                                 df$meds_insulin_sr_ni.3 == 1,
                               1, 0)
table(df$meds_insulin_sr_ni, useNA="always")
df$meds_metformin_sr_ni = ifelse(df$meds_metformin_sr_ni.0 == 1 |
                                   df$meds_metformin_sr_ni.1 == 1 |
                                   df$meds_metformin_sr_ni.2 == 1 |
                                   df$meds_metformin_sr_ni.3 == 1,
                                 1, 0)
table(df$meds_metformin_sr_ni, useNA="always")
df$meds_nonmet_oad_sr_ni = ifelse(df$meds_nonmet_oad_sr_ni.0 == 1 |
                                    df$meds_nonmet_oad_sr_ni.1 == 1 |
                                    df$meds_nonmet_oad_sr_ni.2 == 1 |
                                    df$meds_nonmet_oad_sr_ni.3 == 1,
                                  1, 0)
table(df$meds_nonmet_oad_sr_ni, useNA="always")

# Combination medication variables ----------------------------------------
### remove additions from HSD code and put in ==1 | 

df$meds_any_sr_ni <- 0 #this single variable captures all 3 categories of DM medications
df$meds_any_sr_ni[df$meds_insulin_sr_ni ==1 |  df$meds_metformin_sr_ni ==1 |  df$meds_nonmet_oad_sr_ni ==1 ] <- 1
table(df$meds_any_sr_ni, useNA="always")
table(df$meds_metformin_sr_ni, useNA="always")
table(df$meds_nonmet_oad_sr_ni, useNA="always")	
df$meds_any_sr_ni_ts <- 0 #this single variable captures all 3 categories of DM medications AND the insulin TS
df$meds_any_sr_ni_ts[df$meds_insulin_sr_ni ==1 |  df$meds_metformin_sr_ni ==1 |  df$meds_nonmet_oad_sr_ni ==1 |  df$dm_insulin_sr_ts ==1 ] <- 1
table(df$dm_insulin_sr_ts, useNA="always")

table(df$meds_any_sr_ni, useNA="always")	
table(df$meds_any_sr_ni_ts, useNA="always")


# 1.1  ---------------------------------------------------------------

#1.1 ## 
df$dm_unlikely <- 1
## any SR diabetes or T1DM or T2DM from NI -> they have 0 for diabetes unlikely
df$dm_unlikely[df$dm_any_sr_ni == 1 | df$dm_t1dm_sr_ni == 1 | df$dm_t2dm_sr_ni ==1 ] <- 0
## any on diabetes medications -> they have 0 for diabetes unlikely
df$dm_unlikely[df$meds_any_sr_ni == 1 | df$dm_insulin_sr_ts == 1] <- 0
## any female with gestational diabetes -> they have 0 for diabetes unlikely
df$dm_unlikely[df$`31-0.0`==0 & df$dm_gdm_sr_ni==1 | df$dm_gdm_sr_ts==1] <- 0
## create a field if diabetes is NOT unlikely -> pass on step 1.2
df$tostep1.2 <- ifelse(df$dm_unlikely == 1, 0, 1)
table(df$tostep1.2, useNA = "always")


# 1.2  --------------------------------------------------------------------

#1.2
df$possible_gdm <- 0
## if they are NOT unlikely for diabetes, and they have gestational diabets (TS) and no medication and no insulin and no T1D or T2D diagnosis -> then they have possible gestational diabetes only
df$possible_gdm[df$tostep1.2==1 & df$dm_gdm_sr_ts == 1 & df$meds_any_sr_ni == 0 & df$dm_insulin_sr_ts == 0 & df$dm_t1dm_sr_ni == 0 & df$dm_t2dm_sr_ni == 0] <- 1
## if they are NOTunlikely for diabetes, and they have gestational diabetes (NI) and diagnosis of any diabetes is < 50, and they have no medication and no insulin and and no T1D or T2D diagnosis -> then they have possible gestational diabetes only
df$possible_gdm[df$tostep1.2==1 & df$dm_gdm_sr_ni == 1 & df$age_dm_gdm_dx_sr < 50 & df$meds_any_sr_ni == 0 & df$dm_insulin_sr_ts == 0 & df$dm_t1dm_sr_ni == 0 & df$dm_t2dm_sr_ni == 0] <- 1
df$tostep1.3 <- 0
## create a field if diabetes is NOT unlikely and possible gestational diabetes is 0 -> pass on to step 1.3
## 1.3 = likely diabetes, NOT gestational 
df$tostep1.3[df$tostep1.2 == 1 & df$possible_gdm == 0] <- 1
table(df$tostep1.3, useNA = "always")



# 1.3 ---------------------------------------------------------------------

df$possible_t2dm_a <- 0
## if they might have diabetes and might not have gestational diabets, AND non-metformin OAD
df$possible_t2dm_a[df$tostep1.3==1 & df$meds_nonmet_oad_sr_ni == 1] <- 1
df$tostep1.4 <- 0
## 1.4 = likely diabtes, not gestational, and NOTTT non-metformin OAD 
df$tostep1.4[df$tostep1.3 == 1 & df$possible_t2dm_a == 0] <- 1
table(df$tostep1.4, useNA = "always")

#### age of diagnosis is between 0 and 36. This equals 1. 
df$allancestries_36 <- ifelse(df$agedm_ts_or_ni > 0 & df$agedm_ts_or_ni < 37, 1, 0)


# 1.4 ---------------------------------------------------------------------
df$tostep1.5 <- 0
#### likely diabetes, not gestational, NOT non-metformim OAD, EUR and diagnosed 0-36
df$tostep1.5[df$tostep1.4==1 & df$allancestries_36 == 1] <- 1
df$possible_t2dm_b <- 0
## likely diabetes, not gestational, NOT non-metformin OAD, EUR, diagnosed ABOVE 36 
df$possible_t2dm_b[df$tostep1.4 ==1 & df$tostep1.5 == 0] <- 1
table(df$tostep1.5, useNA = "always")


# 1.5 ---------------------------------------------------------------------

df$possible_t1dm_temp <- 0
#### likely diabetes, young dx, and insulin and/or insulin within 1year of diagnosis
df$possible_t1dm_temp[df$tostep1.5 ==1 & df$dm_insulin_sr_ts == 1 & df$meds_insulin_sr_ni ==1] <-1
df$possible_t1dm_temp[df$tostep1.5 ==1 & (df$`2986-0.0` == 1 | df$`2986-1.0` == 1 | df$`2986-2.0` == 1)] <- 1
## liekly diabetes, young dx, self reported T1D
df$possible_t1dm_temp[df$tostep1.5 ==1 & df$dm_t1dm_sr_ni ==1] <- 1
table(df$possible_t1dm_temp, useNA="always")
df$possible_t2dm_c <- 0
#### likely diabetes, young dx ages, and NOT self-reported T1D
df$possible_t2dm_c[df$tostep1.5 ==1 & df$possible_t1dm_temp == 0] <- 1

df$possible_t2dm_temp <- 0
#### 
df$possible_t2dm_temp[df$possible_t2dm_a == 1 | df$possible_t2dm_b == 1 | df$possible_t2dm_c == 1] <- 1
table(df$possible_t2dm_temp, useNA="always")


# 2.1 ---------------------------------------------------------------------

df$probable_t1dm_a <- 0
####
df$probable_t1dm_a[df$possible_t1dm_temp == 1 & df$dm_t1dm_sr_ni == 1] <- 1
df$tostep2.2 <- 0
df$tostep2.2[df$possible_t1dm_temp ==1 & df$probable_t1dm_a == 0] <- 1

table(df$possible_t1dm_temp, useNA="always")
table(df$probable_t1dm_a, useNA="always")
table(df$tostep2.2, useNA="always")

# 2.2 --------------------------------------------------------------------

df$probable_t1dm_b <- 0
#df$probable_t1dm_b[df$tostep2.2 ==1 & df$dm_insulin_sr_ts == 1 & df$meds_insulin_sr_ni ==1] <- 1
df$probable_t1dm_b[df$tostep2.2 ==1 & df$dm_insulin_sr_ts == 1 & df$`2986-0.0` == 1] <- 1
df$probable_t1dm_b[df$tostep2.2 ==1 & df$meds_insulin_sr_ni == 1 & df$`2986-0.0` == 1] <- 1

df$probable_t1dm <- 0
df$probable_t1dm[df$probable_t1dm_a == 1 | df$probable_t1dm_b ==1] <- 1

df$possible_t1dm <- 0
df$possible_t1dm[df$tostep2.2 ==1 & df$probable_t1dm_b == 0] <- 1
table(df$possible_t1dm, useNA="always")
#### probable_t1dm_a = possible_t1dm_temp = likely diabetes (notgbm) - young dx ages - and NOT on non-metformin OAD medication, and insulin and/or insulin within 1year of diagnosis, and self reported T1D
#### 2.2 = likely diabetes (notgbm), young dx ages, NOT on non-metformin OAD medication, and insulin and/or insulin within 1year of diagnosis, and self reported T1D
#### probable_t1dm_b = 2.2 + insulin
#### possible_t1dm = a or b:


# 3.1 ---------------------------------------------------------------------

df$tostep3.2 <- 0
df$tostep3.2[df$possible_t2dm_temp == 1 & df$meds_metformin_sr_ni == 1 & df$meds_insulin_sr_ni ==0 & df$meds_nonmet_oad_sr_ni ==0 & df$dm_insulin_sr_ts ==0] <- 1
table(df$tostep3.2, useNA = "always")

df$tostep3.3_a <- 0
df$tostep3.3_a[df$possible_t2dm_temp == 1 & df$tostep3.2 == 0] <- 1
table(df$tostep3.3_a, useNA = "always")

# 3.2 ---------------------------------------------------------------------

df$dm_unlikely_b <- 0
df$dm_unlikely_b[df$tostep3.2 ==1 & df$dm_any_sr_ni == 0 & df$dm_gdm_sr_ni == 0 & df$dm_t1dm_sr_ni == 0 & df$dm_t2dm_sr_ni == 0] <- 1
df$dm_unlikely[df$dm_unlikely_b==1] <-1
table(df$dm_unlikely, useNA = "ifany")

df$tostep3.3_b <- 0
df$tostep3.3_b[df$tostep3.2 ==1 & df$dm_unlikely_b == 0] <- 1

df$tostep3.3 <- 0
df$tostep3.3[df$tostep3.3_a == 1 | df$tostep3.3_b ==1] <- 1
table(df$tostep3.3, useNA = "ifany")



# 3.3 ---------------------------------------------------------------------

df$probable_t2dm_a <- 0
df$probable_t2dm_a[df$tostep3.3 == 1 & df$meds_nonmet_oad_sr_ni ==1] <- 1
df$tostep3.4 <- 0
df$tostep3.4[df$tostep3.3 == 1 & df$probable_t2dm_a ==0] <- 1
table(df$tostep3.4, useNA = "ifany")


# 3.4 ---------------------------------------------------------------------

df$probable_t2dm_b <- 0
df$probable_t2dm_b[df$tostep3.4 == 1 & df$dm_insulin_sr_ts == 0 & df$meds_insulin_sr_ni == 0] <- 1

df$probable_t2dm <- 0
df$probable_t2dm[df$probable_t2dm_a == 1 | df$probable_t2dm_b ==1] <- 1
table(df$probable_t2dm, useNA = "ifany")

df$tostep3.5 <- 0
df$tostep3.5[df$tostep3.4 == 1 & df$probable_t2dm == 0] <- 1
table(df$tostep3.5, useNA = "ifany")


# 3.5 ---------------------------------------------------------------------

df$possible_t2dm <- 0
df$possible_t2dm[df$tostep3.5 ==1 & df$dm_t1dm_sr_ni == 0] <- 1
table(df$possible_t2dm, useNA = "ifany")

df$probable_t1dm_c <- 0
df$probable_t1dm_c[df$tostep3.5 ==1 & df$dm_t1dm_sr_ni == 1] <- 1
df$probable_t1dm[df$probable_t1dm_c ==1] <- 1
table(df$probable_t1dm, useNA = "ifany")


# Flag age discrepancy ----------------------------------------------------

### this is anyone who has diagnosis within each category more than 10 of
### gdm, any, t1dm, t2dm, and 2976

df$age_discrepancy_all = NULL
df$age_discrepancy_all = ifelse(
  (
    abs(df$age_dm_any_dx_sr.0 - df$age_dm_any_dx_sr.1) > 10 |
      abs(df$age_dm_any_dx_sr.0 - df$age_dm_any_dx_sr.2) > 10 |
      abs(df$age_dm_any_dx_sr.0 - df$age_dm_any_dx_sr.3) > 10 |
      abs(df$age_dm_any_dx_sr.1 - df$age_dm_any_dx_sr.2) > 10 |
      abs(df$age_dm_any_dx_sr.1 - df$age_dm_any_dx_sr.3) > 10 |
      abs(df$age_dm_any_dx_sr.2 - df$age_dm_any_dx_sr.3) > 10 |
      abs(df$age_dm_gdm_dx_sr.0 - df$age_dm_gdm_dx_sr.1) > 10 |
      abs(df$age_dm_gdm_dx_sr.0 - df$age_dm_gdm_dx_sr.2) > 10 |
      abs(df$age_dm_gdm_dx_sr.0 - df$age_dm_gdm_dx_sr.3) > 10 |
      abs(df$age_dm_gdm_dx_sr.1 - df$age_dm_gdm_dx_sr.2) > 10 |
      abs(df$age_dm_gdm_dx_sr.1 - df$age_dm_gdm_dx_sr.3) > 10 |
      abs(df$age_dm_gdm_dx_sr.2 - df$age_dm_gdm_dx_sr.3) > 10 |
      abs(df$age_dm_t2dm_dx_sr.0 - df$age_dm_t2dm_dx_sr.1) > 10 |
      abs(df$age_dm_t2dm_dx_sr.0 - df$age_dm_t2dm_dx_sr.2) > 10 |
      abs(df$age_dm_t2dm_dx_sr.0 - df$age_dm_t2dm_dx_sr.3) > 10 |
      abs(df$age_dm_t2dm_dx_sr.1 - df$age_dm_t2dm_dx_sr.2) > 10 |
      abs(df$age_dm_t2dm_dx_sr.1 - df$age_dm_t2dm_dx_sr.3) > 10 |
      abs(df$age_dm_t2dm_dx_sr.2 - df$age_dm_t2dm_dx_sr.3) > 10 |
      abs(df$age_dm_t1dm_dx_sr.0 - df$age_dm_t1dm_dx_sr.1) > 10 |
      abs(df$age_dm_t1dm_dx_sr.0 - df$age_dm_t1dm_dx_sr.3) > 10 |
      abs(df$age_dm_t1dm_dx_sr.0 - df$age_dm_t1dm_dx_sr.2) > 10 |
      abs(df$age_dm_t1dm_dx_sr.1 - df$age_dm_t1dm_dx_sr.2) > 10 |
      abs(df$age_dm_t1dm_dx_sr.1 - df$age_dm_t1dm_dx_sr.3) > 10 |
      abs(df$age_dm_t1dm_dx_sr.2 - df$age_dm_t1dm_dx_sr.3) > 10 |
      abs(df$`2976-0.0` - df$`2976-1.0`) > 10 |
      abs(df$`2976-0.0` - df$`2976-2.0`) > 10 |
      abs(df$`2976-0.0` - df$`2976-3.0`) > 10 |
      abs(df$`2976-1.0` - df$`2976-2.0`) > 10 |
      abs(df$`2976-1.0` - df$`2976-3.0`) > 10	|
      abs(df$`2976-2.0` - df$`2976-3.0`) > 10	)		
  , 1, NA)
df$age_discrepancy_all[is.na(df$age_discrepancy_all)] <- 0
table(df$age_discrepancy_all, useNA="always")


# Finishing up ------------------------------------------------------------------
library(data.table)
old = c(
  "dm_unlikely",
  "possible_t1dm",
  "probable_t1dm",
  "possible_t2dm",
  "probable_t2dm",
  "possible_gdm",
  "agedm_ts_or_ni",
  "meds_insulin_sr_ni",
  "meds_metformin_sr_ni",
  "meds_nonmet_oad_sr_ni",
  "meds_any_sr_ni",
  "meds_any_sr_ni_ts",
  "dm_insulin_sr_ts"
)
all(c(old) %in% colnames(df))
setnames(
  df,
  old = c(
    "dm_unlikely",
    "possible_t1dm",
    "probable_t1dm",
    "possible_t2dm",
    "probable_t2dm",
    "possible_gdm",
    "agedm_ts_or_ni",
    "meds_insulin_sr_ni",
    "meds_metformin_sr_ni",
    "meds_nonmet_oad_sr_ni",
    "meds_any_sr_ni",
    "meds_any_sr_ni_ts",
    "dm_insulin_sr_ts"
  ),
  new = c(
    "dm_unlikely_all",
    "possible_t1dm_all",
    "probable_t1dm_all",
    "possible_t2dm_all",
    "probable_t2dm_all",
    "possible_gdm_all",
    "agedm_ts_or_ni_all",
    "meds_insulin_sr_ni_all",
    "meds_metformin_sr_ni_all",
    "meds_nonmet_oad_sr_ni_all",
    "meds_any_sr_ni_all",
    "meds_any_sr_ni_ts_all",
    "dm_insulin_sr_ts_all"
  ),
  skip_absent = TRUE
)

UKB_diabetes <- subset(
  df,
  select = c(
    "eid",
    "2443-0.0",
    "2443-1.0",
    "2443-2.0",
    "2443-3.0",
    "age_discrepancy_all",
    "agedm_ts_or_ni_all",
    "dm_insulin_sr_ts_all",
    "meds_insulin_sr_ni_all",
    "meds_metformin_sr_ni_all",
    "meds_nonmet_oad_sr_ni_all",
    "meds_any_sr_ni_all",
    "meds_any_sr_ni_ts_all",
    "dm_unlikely_all",
    "possible_t1dm_all",
    "probable_t1dm_all",
    "possible_t2dm_all",
    "probable_t2dm_all",
    "possible_gdm_all"
  )
) %>% rename(id = eid)

UKB_diabetes <- UKB_diabetes %>%
  left_join(a1c %>% select(id, hba1c.30750.NGSP.max), by = "id") %>%
  left_join(
    icd_codes %>% select(
      id,
      E10_main,
      E10_secondary,
      icd9_main,
      icd9_secondary,
      icd10_main,
      icd10_secondary
    ),
    by = "id"
  ) %>%
  left_join(diabetes_df %>% select(id = eid, father_dm, mother_dm, sibling_dm),
            by = "id")

# T1D Definition ----------------------------------------------------------
# Case: Possible t1dm OR Probable t1dm

UKB_diabetes <- UKB_diabetes %>%
  mutate(t1d_case = ifelse(possible_t1dm_all == 1 |
                             probable_t1dm_all == 1, 1, 0))
table(UKB_diabetes$t1d_case, useNA = "always")

# T2D Definition ----------------------------------------------------------

# Possible t2dm OR Probable t2dm OR [(max hba1c >= 6.5 & possible t1dm!=1 & probable t1dm!=1 & icd10 E10!=1)]

UKB_diabetes <- UKB_diabetes %>% mutate(
  t2d_case = case_when(
    possible_t2dm_all == 1 ~ 1,
    probable_t2dm_all == 1 ~ 1,
    hba1c.30750.NGSP.max >= 6.5 &
      possible_t1dm_all != 1 &
      probable_t1dm_all != 1 & E10_main != 1 & E10_secondary != 1 ~ 1,
    TRUE ~ 0
  )
)
table(UKB_diabetes$t2d_case, useNA = "always")
table(UKB_diabetes$t2d_case, UKB_diabetes$t1d_case, useNA = "always")


# Diabetes Control --------------------------------------------------------

# Control:
# Dm unlikely
# AND does not have:
# ICD9 or ICD10 codes for any diabetes
# DM medication use: Insulin, metformin, non-metformin OAD (meds_any_sr_ni_ts_all)
# 2443 self-report diabetes (2443-0.0, 2443-1.0, 2443-2.0, 2443-3.0)
# Family history of diabetes (20107, 20110, 20111)
# & max hba1c is < 5.7 (hba1c.30750.NGSP.max)

UKB_diabetes <- UKB_diabetes %>% mutate(
  diabetes_control = case_when(
    dm_unlikely_all == 1 &
      meds_any_sr_ni_ts_all != 1 &
      sibling_dm != 1 &
      mother_dm != 1 &
      father_dm != 1 &
      `2443-0.0` %in% c(-3,-1, 0, NA) &
      `2443-1.0` %in% c(-3,-1, 0, NA) &
      `2443-2.0` %in% c(-3,-1, 0, NA) &
      `2443-3.0` %in% c(-3,-1, 0, NA) &
      hba1c.30750.NGSP.max < 5.7 &
      icd10_main != 1 &
      icd10_secondary != 1 &
      icd9_main != 1 &
      icd9_secondary != 1 ~ 1,
    TRUE ~ 0
  )
)

table(UKB_diabetes$diabetes_control, useNA = "always")
table(UKB_diabetes$diabetes_control, UKB_diabetes$t2d_case, useNA = "always")
table(UKB_diabetes$diabetes_control, UKB_diabetes$t1d_case, useNA = "always")

# Write file --------------------------------------------------------------


write_csv(UKB_diabetes, filename, na = "")
