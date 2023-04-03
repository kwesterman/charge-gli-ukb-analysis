########## Load libraries------------------------------------------------------------------------------------------------------
library(tidyverse)
library(data.table)


########List of fields----------------------------------------------------------------------------------------------------------
all_fields <- c(2090,2100,20002,20003,20126,4598,4609,4620,4631,5375,5386,41202,41204,20446,20441,20436,20439,20440,20449,20536,20532,20435,20450,20437,4642,4653,6156,5663,5674,20544,20514,20510,20534,20519,20511,20507,20508,20518,20513,20400)


mhq_fields <- c(20500, 20499, 20544, 20514, 20510, 20534, 20519, 20511, 20507, 20508, 20518, 20513, 20446, 20441, 20447, 20436, 20439, 20449, 20536, 20532, 20533, 
                20535, 20534, 20435, 20450, 20437, 20438, 20440, 20442, 20433, 20445, 20434, 20448, 20546, 20547, 20501, 20502, 20548, 20492, 20493, 20506, 20509, 20520, 20515, 20516, 20505, 20512, 20412, 20420, 20425, 20425,
                20542, 20538, 20543, 20541, 20540, 20539, 20537, 20426, 20423, 20429, 20419, 20422, 20417, 20427, 20428, 20549, 20550, 20418, 20401, 20406, 20415, 20404, 20503, 20551, 20504, 20456, 20457, 20431, 20552, 20432,
                20414, 20403, 20416, 20413, 20407, 20412, 20409, 20408, 20411, 20405, 20454, 20455, 20471, 20473, 20463, 20465, 20474, 20476, 20468, 20470, 20467, 20461, 20462, 20477, 20466, 20489, 20488, 20487, 20490, 20491, 
                20522, 20523, 20521, 20524, 20525, 20531, 20529, 20526, 20530, 20528, 20527, 20497, 20498, 20495, 20496, 20494, 20479, 20485, 20486, 20480, 20482, 20481, 20553, 20483, 20484, 20554, 20458, 20459, 20460)

phq_fields <-c('20514', '20510', '20517', '20519', '20511', '20507', '20508', '20518', '20513')
all_fields_str <- as.character(all_fields)
mhq_fields_str <- as.character(mhq_fields)

########### Files being used-------------------------------------------------------------------------------------------------------

mental_healthfile1 <-"../data/raw/ukb671261.tab.gz"
mental_healthfile2 <-"../data/raw/ukb47880.csv"

######### Load the files and merge--------------------------------------------------------------------------------------------------
df1 <- fread(mental_healthfile1,data.table=FALSE)

df2_cols <- fread(mental_healthfile2, nrows = 1,data.table=FALSE)
col <- names(df2_cols)

substring2 <- grep(paste(all_fields_str, collapse = "|"), col, value = TRUE)
substring2 <- c(substring2, "eid")
df2 <- fread(mental_healthfile2,select =substring2, data.table=FALSE)
dim(df2)

names(df2) <- paste0("f.", gsub("-", ".", names(df2), fixed = TRUE))

df1 <- df1[order(names(df1))] # sort by column names
df2 <- df2[order(names(df2))]
df <- merge(df1, df2, by = "f.eid")
dim(df)



########### Create has_MHQ column---------------------------------------------------------------------------------------------------

#df[grep("20400", names(df))] # only has one instance
##fun <- function(df){
df$has_mhq <- ifelse(is.na(df[['f.20400.0.0']]), 0, 1)
table(df$has_mhq)

############### Check if any MHQ values are missing----------------------------------------------------------------------------------



###################### define Psychosis------------------------------------------------------------------------------------------------

antipsychotic_meds = list(1140868170, 1140928916, 1141152848, 1140867444, 1140879658, 1140868120, 1141153490, 1140867304, 1141152860, 1140867168, 
1141195974, 1140867244, 1140867152, 1140909800, 1140867420, 1140879746, 1141177762, 1140867456, 1140867952, 1140867150, 1141167976, 1140882100, 1140867342, 
1140863416, 1141202024, 1140882098, 1140867184, 1140867092, 1140882320, 1140910358, 1140867208, 1140909802, 1140867134, 
1140867306, 1140867210, 1140867398, 1140867078, 1140867218, 1141201792, 1141200458, 1140867136, 1140879750, 1140867180, 1140867546, 1140928260, 1140927956)


delusional_disorders_icd10 = list("F20", "F200", "F201", "F202", "F203", "F204", "F205", "F206", "F208", "F209", "F21", "F22", "F220", 
"F228", "F229", "F23", "F230", "F231", "F232", "F233", "F238", "F239", "F24", "F25", "F250", "F251", "F252", "F258", "F259", "F28", "F29")
mood_disorders_no_depression_icd10 = list("F30", "F300", "F301", "F302", "F308", "F309", "F31", "F310", "F311", "F312", "F313", "F314", "F315", "F316", "F317", "F318", "F319", "F34",
 "F340", "F341", "F348", "F349", "F38", "F380", "F381", "F388", "F39")

psychosis_icd10 = append(delusional_disorders_icd10, mood_disorders_no_depression_icd10)

mhq_psychosis_codes = list(2,3,10)


df<-df %>% 
  mutate("6150.0" = rowSums(across(starts_with("f.6156.0.")), na.rm=TRUE))

df<-df %>% 
  mutate("6150.1" = rowSums(across(starts_with("f.6156.1.")), na.rm=TRUE))

df<-df %>% 
  mutate("6150.2" = rowSums(across(starts_with("f.6156.2.")), na.rm=TRUE))

df<-df %>% 
  mutate("6150.3" = rowSums(across(starts_with("f.6156.3.")), na.rm=TRUE))
  

table(df[["6150.0"]])


for (i in 0:3)
{
  # create flag columns with default value 0
  #  df[[paste0('x', toString(1))]] = 100
  df[[paste0('sr_psychosis_',toString(i))]]=0
  df[[paste0('psychosis_meds_',toString(i))]]=0
  df[[paste0('hospital_psychosis_',toString(i))]]=0
  df[[paste0('mhq_psychosis_',toString(i))]]=0
  df[[paste0('hyper_manic_symptoms',toString(i))]]=0
  
  # based on the condition below, change flag values to 1 if satisfied
  # 1. self reported psychosis  
  df[(df[[paste0('f.20002.',toString(i),'.0')]] %in% 1289)| (df[[paste0('f.20002.',toString(i),'.0')]] %in%1291),paste0('sr_psychosis_',toString(i))] <-1
  # 2. Antipsychotic Usage
  df[df[[paste0('f.20003.',toString(i),'.0')]] %in% antipsychotic_meds, paste0('psychosis_meds_',toString(i))] <-1
  # 3. Hospital ICD 10 Psychosis
  df[(df[[paste0('f.41202.',toString(i),'.0')]] %in% psychosis_icd10)| (df[[paste0('f.41204.',toString(i),'.0')]] %in% psychosis_icd10),paste0('hospital_psychosis_',toString(i))] <-1
  # 4. Hypermanic Symptoms
  df[(df[[paste0('6150.',toString(i))]]==15)| (df[[paste0('6150.',toString(i))]]>30),paste0('hyper_manic_symptoms_',toString(i))] <-1
  # 5. MHQ Psychosis
  df[df[[paste0('f.20544.',toString(i),'.0')]] %in% mhq_psychosis_codes, paste0('mhq_psychosis_',toString(i))] <-1
  
}  

for (i in 0:3)
{
   
   # create MD variables and assign default value 0
   df[[paste0('bipolar_type1_',toString(i))]] =0
   df[[paste0('bipolar_type2_',toString(i))]] =0
   df[[paste0('is_bipolar_',toString(i))]] =0
      
   
   df[((df[[paste0('4642', toString(i), '.0')]]==1)|(df[[paste0('4653', toString(i), ".0")]]==1))&(df[[paste0('hyper_manic_symptoms', toString(i))]]==1) &(df[[paste0('5663', toString(i), ".0")]]==1)&(df[[paste0('5674', toString(i), ".0")]]==1), paste0('bipolar_type1_', toString(i))]=1
   
   
   df[((df[[paste0('4642', toString(i), '.0')]]==1)|(df[[paste0('4653', toString(i), ".0")]]==1))&(df[[paste0('hyper_manic_symptoms', toString(i))]]==1) &(df[[paste0('5663', toString(i), ".0")]]==1), paste0('bipolar_type2_', toString(i))]=1
   
   df[(df[[paste0('bipolar_type1_', toString(i), '.0')]]==1)|(df[[paste0('bipolar_type2_', toString(i), ".0")]]==1), paste0('is_bipolar_', toString(i))]=1
   
 }

for (i in 0:3)
{

df[[paste0('has_psychosis',toString(i))]]=0  

  ### for participants with MHQ:
  #### (has_mhq==1) &[(sr_psychosis ==1) | (psychosis_meds ==1 ) | (mhq_psychosis == 1) | (hospital_psychosis ==1)  | (is_bipolar=1)]
  df[(df$has_mhq ==1)&((df[[paste0('sr_psychosis_', toString(i))]]==1)|(df[[paste0('psychosis_meds_', toString(i))]]==1)|(df[[paste0('mhq_psychosis_', toString(i))]]==1) |(df[[paste0('hospital_psychosis_', toString(i))]]==1)|(df[[paste0('is_bipolar_', toString(i))]]==1) | (df[[paste0('20126', toString(i), ".0")]]==1)), paste0('has_psychosis', toString(i))]=1
 
  ### for participants without MHQ:
  #### [(has_mhq ==0) &(sr_psychosis ==1) | (psychosis_meds ==1 ) | (hospital_psychosis ==1)  | (is_bipolar=1)]
  df[(df$has_mhq ==0)&((df[[paste0('sr_psychosis_', toString(i))]]==1)|(df[[paste0('psychosis_meds_', toString(i))]]==1) |(df[[paste0('hospital_psychosis_', toString(i))]]==1)|(df[[paste0('is_bipolar_', toString(i))]]==1)), paste0('has_psychosis', toString(i))]=1
}

df <- df %>% mutate("has_psychosis_ever" = case_when((rowSums(across(starts_with("has_psychosis"))))>=1~1, (rowSums(across(starts_with("has_psychosis"))))<1~0))
################################################################################################################################################################

################# PHQ score----------------------------------------------------------------------------------------------------------------

phq_subs <- grep(paste(phq_fields, collapse = "|"), colnames(df), value = TRUE)
#df_phq <-df[,select.phq_subs]

df[,c(phq_subs)][is.na(df[,c(phq_subs)])] <- 0
df[,c(phq_subs)][df[,c(phq_subs)] == -818]<-0

df[,c(phq_subs)]<-df[,c(phq_subs)] %>%
  mutate(across(1:last_col(), ~ifelse(.>0, . - 1, .)))
  
df$PHQ_score <- rowSums(df[,c(phq_subs)])





################################# FOR MHQ PARTICIPANTS---------------------------------------------------------------------------------------

for (i in 0:3)
{
  # create flag columns with default value 0
  #  df[[paste0('x', toString(1))]] = 100
  df[[paste0('MHQ_depression_core_',toString(i))]]=0
  df[[paste0('MHQ_threshold_score_',toString(i))]]=0
  df[[paste0('flag_20449_more_tired_',toString(i))]]=0
  df[[paste0('flag_20536_weight_fluctuation_',toString(i))]]=0
  df[[paste0('flag_20532_change_in_sleep_',toString(i))]]=0
  df[[paste0('flag_20435_trouble_concentrating_',toString(i))]]=0
  df[[paste0('flag_20450_self_deprecating_feelings_',toString(i))]]=0
  df[[paste0('flag_20437_thoughts_of_death_',toString(i))]]=0
  #df[[paste0('MHQ_symptoms_for_atleast_5days_',toString(i))]]=0
  # based on the condition below, change flag values to 1 if satisfied
  # 1. MHQ depression core
  df[(df[[paste0('f.20446.',toString(i),'.0')]]%in%1)|(df[[paste0('f.20441.',toString(i),'.0')]]%in%1), paste0('MHQ_depression_core_',toString(i))] <-1
  # 2. MHQ threshold score
  df[(!is.na(df[[paste0('f.20436.',toString(i),'.0')]]))&(df[[paste0('f.20436.',toString(i),'.0')]]>1)&(df[[paste0('f.20439.',toString(i),'.0')]]>1)&(df[[paste0('f.20440.',toString(i),'.0')]]>1),paste0('MHQ_threshold_score_',toString(i))] <- 1
  # 3. More tired
  df[df[[paste0('f.20449.',toString(i),'.0')]]%in%1,paste0('flag_20449_more_tired_',toString(i))]= 1 
  # 4. Fluctuation in weight
  df[(!is.na(df[[paste0('f.20536.',toString(i),'.0')]]))&df[[paste0('f.20536.',toString(i),'.0')]]>0,paste0('flag_20536_weight_fluctuation_',toString(i))]= 1 
  # 5. chnage in sleep patterns
  df[df[[paste0('f.20532.',toString(i),'.0')]]%in%1,paste0('flag_20532_change_in_sleep_',toString(i))]= 1 
  # 6. trouble concentrating
  df[df[[paste0('f.20435.',toString(i),'.0')]]%in%1,paste0('flag_20435_trouble_concentrating_',toString(i))]= 1 
  # 7. self deprecation 
  df[df[[paste0('f.20450.',toString(i),'.0')]]%in%1,paste0('flag_20450_self_deprecating_feelings_',toString(i))]= 1 
  # 8. thoughts of death
  df[df[[paste0('f.20437.',toString(i),'.0')]]%in%1,paste0('flag_20437_thoughts_of_death_',toString(i))]= 1 
   
}  

for(i in 0:3)
{
  df<-df %>%
  mutate(!!paste0('MHQ_symptoms_for_atleast_5days_', toString(i)) := rowSums(across(c(paste0('flag_20449_more_tired_',toString(i)), 
                                                                                  paste0('flag_20536_weight_fluctuation_',toString(i)),
																				  paste0('flag_20532_change_in_sleep_',toString(i)),
																				  paste0('flag_20435_trouble_concentrating_',toString(i)),
																				  paste0('flag_20450_self_deprecating_feelings_',toString(i)),
																				  paste0('flag_20437_thoughts_of_death_',toString(i))
																				  ))))
																				  
																				  
  #df<-df %>% mutate(!!paste0('dDEPR_MHQ_life_time_', toString(i))  := (paste0('MHQ_depression_core_',toString(i))==1)&(paste0('MHQ_threshold_score_',toString(i))==1)&(paste0('MHQ_symptoms_for_atleast_5days_', toString(i))>=5))
}  


############################ For non-MHQ participants---------------------------------------------------------------------------------------------------

antidepressant_meds = list(1140879616, 1140921600, 1140879540, 1140867878, 1140916282, 1140909806, 1140867888, 1141152732, 1141180212, 1140879634, 1140867876, 
1140882236, 1141190158, 1141200564, 1140867726,  1140879620, 1140867818, 1140879630, 1140879628, 1141151946, 1140867948, 1140867624,1140867756, 1140867884, 
1141151978, 1141152736, 1141201834, 1140867690, 1140867640, 1140867920, 1140867850, 1140879544, 1141200570, 1140867934, 1140867758, 1140867914, 1140867820, 
1141151982, 1140882244, 1140879556, 1140867852, 1140867860, 1140917460, 1140867938, 1140867856, 1140867922, 1140910820, 1140882312, 1140867944, 1140867784, 1140867812, 1140867668)


depressive_episodes_icd10 = list("F32", "F320", "F321", "F322", "F323", "F328", "F329")
recurrent_depressive_disorder_icd10 = list("F33", "F330", "F331", "F332", "F333", "F334", "F338", "F339")

depression_icd10 = append(depressive_episodes_icd10, recurrent_depressive_disorder_icd10)


for (i in 0:3)
{
  # create flag columns with default value 0
  #  df[[paste0('x', toString(1))]] = 100
  df[[paste0('help_seeking_',toString(i))]]=0
  df[[paste0('sr_depression_',toString(i))]]=0
  df[[paste0('sr_antidepressant_usage_',toString(i))]]=0
  df[[paste0('depression_hospital_icd10_',toString(i))]]=0
  
  
  # 1. help seeking
  df[(df[[paste0('f.2090.',toString(i),'.0')]]%in%1)|(df[[paste0('f.2100.',toString(i),'.0')]]%in%1), paste0('help_seeking_',toString(i))] <-1
  # 2. self reported depression
  df[(df[[paste0('f.20002.',toString(i),'.0')]] %in% 1286), paste0('sr_depression_',toString(i))] <- 1
  # 3. Anti depressant usage - Self reported
  df[df[[paste0('f.20003.',toString(i),'.0')]] %in% antidepressant_meds, paste0('sr_antidepressant_usage_',toString(i))] <-1 
  # 4. ICD 10 depression
  df[(df[[paste0('f.41202.',toString(i),'.0')]] %in% depression_icd10)| (df[[paste0('f.41204.',toString(i),'.0')]] %in% depression_icd10),paste0('depression_hospital_icd10_',toString(i))] <-1
     
}  

# Smith definition of depression
df[df[['f.20126.0.0']]%in%5, 'derived_single_md_0'] <-1
df[df[['f.20126.0.0']]%in%4, 'derived_repeat_moderate_md_0'] <-1
df[df[['f.20126.0.0']]%in%3, 'derived_repeat_severe_md_0'] <-1

# create flag variable columns for MD definitions
for (i in 1:3)
{
   # create MD variables and assign default value 0
   df[[paste0('derived_single_md_',toString(i))]] =0
   df[[paste0('derived_repeat_moderate_md_',toString(i))]] =0
   df[[paste0('derived_repeat_severe_md_',toString(i))]] =0
   
   # Single probable MD
   df[((df[[paste0('f.4598.',toString(i),'.0')]]%in%1)& (df[[paste0('f.4609.',toString(i),'.0')]]>=2)& (df[[paste0('f.4620.',toString(i),'.0')]]%in%1) &((df[[paste0('f.2090.',toString(i),'.0')]]%in%1) | (df[[paste0('f.2100.',toString(i),'.0')]]%in%1)))|
   ((df[[paste0('f.4631.',toString(i),'.0')]]%in%1)& (df[[paste0('f.5375.',toString(i),'.0')]]>=2)& (df[[paste0('f.5386.',toString(i),'.0')]]%in%1) &((df[[paste0('f.2090.',toString(i),'.0')]]%in%1) | (df[[paste0('f.2100.',toString(i),'.0')]]%in%1))), 
   paste0('derived_single_md_',toString(i))] <- 1
   
   # Probable Recurrent MD (moderate)
   df[((df[[paste0('f.4598.',toString(i),'.0')]]%in%1)&(df[[paste0('f.4609.',toString(i),'.0')]]>=2)&(df[[paste0('f.4620.',toString(i),'.0')]]%in%1)&(df[[paste0('f.2090.',toString(i),'.0')]]%in%1))|
   ((df[[paste0('f.4631.',toString(i),'.0')]]%in%1)&(df[[paste0('f.5375.',toString(i),'.0')]]>=21)&(df[[paste0('f.5386.',toString(i),'.0')]]%in%1)&(df[[paste0('f.2090.',toString(i),'.0')]]%in%1)), 
   paste0('derived_repeat_moderate_md_', toString(i))]=1
   
   # Probable Recurrent MD (severe): [(1) OR (2)] AND (4) AND (5) AND (7)
   df[((df[[paste0('f.4598.',toString(i),'.0')]]%in%1)& (df[[paste0('f.4609.',toString(i),'.0')]]>=2) &(df[[paste0('f.4620.',toString(i),'.0')]]%in%1) & (df[[paste0('f.2100.',toString(i),'.0')]]%in%1))|
   ((df[[paste0('f.4631.',toString(i),'.0')]]%in%1) & (df[[paste0('f.5375.',toString(i),'.0')]]>=2)& (df[[paste0('f.5386.',toString(i),'.0')]]%in%1) & (df[[paste0('f.2100.',toString(i),'.0')]]%in%1)), 
   paste0('derived_repeat_severe_md_', toString(i))]=1
}

### create a subset with ONLY the the definition scripts and save in saveRDS
subset_cols = c('f.eid', 'has_psychosis0', 'has_mhq', 'PHQ_score', 'MHQ_depression_core_0', 'MHQ_threshold_score_0', 'MHQ_symptoms_for_atleast_5days_0', 'help_seeking_0', 'sr_depression_0', 'sr_antidepressant_usage_0',
'depression_hospital_icd10_0', 'derived_single_md_0', 'derived_repeat_moderate_md_0', 'derived_repeat_severe_md_0')

df_final = df[, subset_cols]

saveRDS(df_final, "depression_intermediate_df.rds")

###################################
### following function filters, residualizes and scales to calcuate the qDEPR value
##################################
calculate_qDEPR <- function(df, y) {
  qDEPR_mask <- (df$has_mhq == 1) & (df$has_psychosis0 == 0)
  lm_str <- paste0(y, " ~ sex * age")
  lm_fit <- lm(as.formula(lm_str), data=df[qDEPR_mask, ], na.action=na.exclude)
  df$qDEPR <- NA
  df$qDEPR[qDEPR_mask] <- as.vector(scale(resid(lm_fit)))
  df
}