library(tidyverse)
library(data.table)

filename <- "ukb_glycemictraits.csv"

input_file1 <- "../data/raw/ukb51002.csv"
input_file2 <- "../data/raw/ukb50493.csv"


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
input_file1 <- fread(input_file1, data.table = FALSE) %>%
  filter(eid %in% valid_ids)
input_file2 <- fread(input_file2, data.table = FALSE) %>%
  filter(eid %in% valid_ids)


# Create dataframe --------------------------------------------------------


glycemic_df <- input_file1 %>%
  select(id = eid,
         hba1c = `30750-0.0`,
         glucose = `30740-0.0`) %>%
  left_join(input_file2 %>%
              select(id = eid,
                     fasting_hours = `74-0.0`),
            by = "id") %>%
  mutate(
    hba1c = hba1c / 10.929 + 2.15,
    glucose = glucose / 0.06,
    fasting_status = ifelse(fasting_hours >= 8, "F", "NF")
  ) %>%
  select(id, hba1c, glucose, fasting_status)


### Winsorization ------------------------------------------------------------
##
### hba1c
##hba1c_upper <-
##  mean(glycemic_df$hba1c, na.rm = T) + (6 * sd(glycemic_df$hba1c, na.rm = T))
##hba1c_lower <-
##  mean(glycemic_df$hba1c, na.rm = T) - (6 * sd(glycemic_df$hba1c, na.rm = T))
##table(glycemic_df$hba1c > hba1c_upper) # 1748 above
##table(glycemic_df$hba1c < hba1c_lower) # none below
##
##table(glycemic_df$fasting_status, useNA = "ifany")
##
### glucose - fasted (F)
##fasted_upper <-
##  mean(glycemic_df %>% filter(fasting_status == "F") %>% pull(glucose),
##       na.rm = T) +
##  (6 * sd(
##    glycemic_df %>% filter(fasting_status == "F") %>% pull(glucose),
##    na.rm = T
##  ))
##fasted_lower <-
##  mean(glycemic_df %>% filter(fasting_status == "F") %>% pull(glucose),
##       na.rm = T) -
##  (6 * sd(
##    glycemic_df %>% filter(fasting_status == "F") %>% pull(glucose),
##    na.rm = T
##  ))
##print(paste(
##  glycemic_df %>% filter(fasting_status == "F" &
##                           glucose > fasted_upper) %>% nrow,
##  "values winsorized fasted upper"
##)) # 89 above
##print(paste(
##  glycemic_df %>% filter(fasting_status == "F" &
##                           glucose < fasted_lower) %>% nrow,
##  "values winsorized fasted lower"
##)) # 0
##
### glucose - non-fasted (NF)
##nonfasted_upper <-
##  mean(glycemic_df %>% filter(fasting_status == "NF") %>% pull(glucose),
##       na.rm = T) +
##  (6 * sd(
##    glycemic_df %>% filter(fasting_status == "NF") %>% pull(glucose),
##    na.rm = T
##  ))
##nonfasted_lower <-
##  mean(glycemic_df %>% filter(fasting_status == "NF") %>% pull(glucose),
##       na.rm = T) -
##  (6 * sd(
##    glycemic_df %>% filter(fasting_status == "NF") %>% pull(glucose),
##    na.rm = T
##  ))
##print(paste(
##  glycemic_df %>% filter(fasting_status == "NF" &
##                           glucose > nonfasted_upper) %>% nrow,
##  "values winsorized non-fasted lower"
##)) # 2157 above
##print(paste(
##  glycemic_df %>% filter(fasting_status == "NF" &
##                           glucose < nonfasted_lower) %>% nrow,
##  "values winsorized non-fasted lower"
##)) # 0
##
##
##glycemic_df <- glycemic_df %>%
##  mutate(
##    hba1c = ifelse(hba1c < hba1c_lower, hba1c_lower, hba1c),
##    hba1c = ifelse(hba1c > hba1c_upper, hba1c_upper, hba1c)
##  )
##
##glycemic_df <- glycemic_df %>%
##  mutate(
##    glucose = ifelse(
##      fasting_status == "F" &
##        glucose < fasted_lower,
##      fasted_lower,
##      glucose
##    ),
##    glucose = ifelse(
##      fasting_status == "F" &
##        glucose > fasted_upper,
##      fasted_upper,
##      glucose
##    ),
##    glucose = ifelse(
##      fasting_status == "NF" &
##        glucose < nonfasted_lower,
##      nonfasted_lower,
##      glucose
##    ),
##    glucose = ifelse(
##      fasting_status == "NF" &
##        glucose > nonfasted_upper,
##      nonfasted_upper,
##      glucose
##    )
##  )


# Write file --------------------------------------------------------------


write_csv(glycemic_df, filename, na = "")
