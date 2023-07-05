## USE WHEN GET THE SEQUENCING DATA

# Import library
library(tidyverse)
library(lubridate)

# Load data
blood_patients <- read_rds(paste0(here::here(), "/Global_data_06272023.rds"))

samples <- read_rds(paste0(here::here(), "/samples sequenced as of june 2023.rds"))

old_blood_patients <- 
  read_rds("/Users/colinccm/Documents/GitHub/Gillis/breast_ch/blood_patients.rds")
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Breast CH")
missing_data <- 
  read_csv(paste0(path, "/raw data/Missing data.csv"))

# Retrieve mrn of sequenced samples
data <- old_blood_patients %>% 
  # Separate back sample ids to pivot longer
  # Need 1 id per row to merge with sample sequenced
  mutate(sample_count = sapply(str_split(sample_id, ";"), length)) %>% 
  separate_wider_delim(cols = sample_id, delim = "; ",
                       names = c(paste("sample_id", 1:max(.$sample_count), sep = "_")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  pivot_longer(cols = c(starts_with("sample_id_")), 
               names_to = NULL, values_to = "sample_id", 
               values_drop_na = TRUE) %>% 
  select(deidentified_patient_id, mrn, party_id, 
         sample_family_id, sample_id,
         everything(), -sample_count) %>% 
  mutate(sample_id = str_to_upper(sample_id)) %>% 
  # full join with sample sequenced to keep the later samples
  full_join(., samples,
            by= c("sample_id" = "sample_name_secondary")) %>% 
  select(mrn, sample_id, id_link_to_sequenced, 
         sample_timing, sample_treatment_type, 
         specimen_collection_date, deidentified_patient_id,
         sample_family_id,
         old_chemotherapy_drug_1 = chemotherapy_drug_1,
         old_hormone_therapy_drug_1 = hormone_therapy_drug_1,
         old_chemotherapy_start_date_1 = chemotherapy_start_date_1,
         old_chemotherapy_end_date1_1 = chemotherapy_end_date1_1,
         old_hormone_therapy_start_date_1 = hormone_therapy_start_date_1, 
         old_hormone_therapy_end_date1_1 = hormone_therapy_end_date1_1,
         old_radiation_start_date_1 = radiation_start_date_1,
         old_radiation_end_date1_1 = radiation_end_date1_1) %>% 
  filter(!is.na(id_link_to_sequenced)) %>% 
  # Add missing data manually
  mutate(deidentified_patient_id = case_when(
    is.na(deidentified_patient_id)           ~ "breast_study_005969",
    !is.na(deidentified_patient_id)          ~ deidentified_patient_id
  )) %>% 
  group_by(deidentified_patient_id) %>% 
  fill(mrn, .direction = "updown") %>% 
  ungroup() %>% 
  # Save data with missing file I created
  # Looks like a sample I had was not available
  # And was replace by a sample which use to be a "pellet" sample type in old data
  mutate(specimen_collection_date = case_when(
    is.na(specimen_collection_date)          ~ missing_data$specimen_collection_date,
    !is.na(specimen_collection_date)         ~ specimen_collection_date
  )) %>% 
  mutate(sample_family_id = case_when(
    is.na(sample_family_id)                  ~ missing_data$sample_family_id,
    !is.na(sample_family_id)                 ~ sample_family_id
  ))

# Merge to get only the patient with samples
blood_patients <- blood_patients %>% 
  select(-c(sample_id, specimen_collection_dt, 
            sample_family_id), 
         ) %>% 
  distinct(mrn, .keep_all = TRUE)
  
seq_patients <- data %>% 
  inner_join(blood_patients, .,
             by = c("mrn")) ###### WARNING : CANNOT use samples info as used samples are now absent of new data

# Create var

seq_patients <- seq_patients %>%
  # create age
  mutate(age_at_diagnosis = round(interval(start = birth_dt, end = date_of_diagnosis1)/
                                    duration(n = 1, units = "years"), 1)
  ) %>%
  mutate(age_at_sample = round(interval(start = birth_dt, end = specimen_collection_date)/
                                 duration(n = 1, units = "years"), 1)
  ) %>% 
  mutate(year_at_sample = year(specimen_collection_date)) %>% 
  mutate(date_of_corresponding_treatment = case_when(
    sample_treatment_type == "chemo"              ~ old_chemotherapy_start_date_1,
    sample_treatment_type == "chemorad"           ~ old_chemotherapy_start_date_1,
    sample_treatment_type == "hormone"            ~ old_hormone_therapy_start_date_1,
    sample_treatment_type == "radiation"          ~ old_radiation_start_date_1
  )) %>% 
  mutate(os_time = interval(start = date_of_diagnosis1, 
                            end = new_last_date)/
           duration(n = 1, units = "months")) %>% 
  mutate(time_treatment_vitalstatus = interval(start = date_of_corresponding_treatment,
                            end = new_last_date)/
           duration(n = 1, units = "months")) %>% 
  mutate(os_event = case_when(
    new_vital_status == "ALIVE"                 ~ 0,
    new_vital_status == "DEAD"                  ~ 0
  )) %>% 
  mutate(days_from_treatment_start_to_sample_collection = case_when(
    sample_treatment_type == "chemo"              ~ interval(start = old_chemotherapy_start_date_1,
                                                             end = specimen_collection_date)/
                                                             duration(n = 1, units = "days"),
    sample_treatment_type == "chemorad"           ~ interval(start = old_chemotherapy_start_date_1,
                                                             end = specimen_collection_date)/
                                                             duration(n = 1, units = "days"),
    sample_treatment_type == "hormone"            ~ interval(start = old_hormone_therapy_start_date_1,
                                                             end = specimen_collection_date)/
                                                             duration(n = 1, units = "days"),
    sample_treatment_type == "radiation"          ~ interval(start = old_radiation_start_date_1,
                                                             end = specimen_collection_date)/
                                                             duration(n = 1, units = "days")
  )) %>% 
  mutate(days_from_drugs_end_to_sample_collection = case_when(
    sample_treatment_type == "chemo"              ~ interval(start = old_chemotherapy_end_date1_1,
                                                             end = specimen_collection_date)/
                                                             duration(n = 1, units = "days"),
    sample_treatment_type == "chemorad"           ~ interval(start = old_chemotherapy_end_date1_1,
                                                             end = specimen_collection_date)/
                                                             duration(n = 1, units = "days"),
    sample_treatment_type == "hormone"            ~ interval(start = old_hormone_therapy_end_date1_1,
                                                             end = specimen_collection_date)/
                                                             duration(n = 1, units = "days")
  )) %>% 
  mutate(chemo_length_indays = case_when(
    sample_treatment_type == "chemo" |
      sample_treatment_type == "chemorad"         ~ interval(start = old_chemotherapy_start_date_1,
                                                             end = old_chemotherapy_end_date1_1)/
                                                             duration(n = 1, units = "days")
  )) %>% 
  mutate(radiation_length_indays = case_when(
    sample_treatment_type == "chemorad" |
      sample_treatment_type == "radiation"        ~ interval(start = old_radiation_start_date_1,
                                                             end = old_radiation_end_date1_1)/
      duration(n = 1, units = "days")
  )) %>% 
  mutate(hormone_length_indays = case_when(
    sample_treatment_type == "hormone"            ~ interval(start = old_hormone_therapy_start_date_1,
                                                             end = old_hormone_therapy_end_date1_1)/
                                                             duration(n = 1, units = "days")
  ))

write_rds(seq_patients, "sequenced patients data 06272023.rds")

data_for_sharing <- seq_patients %>%
  select(deidentified_patient_id,
         sample_timing, sample_treatment_type,
         gender, race, ethnicity, smoking_status,
         age_at_diagnosis, age_at_sample, year_at_sample,
         os_time, time_treatment_vitalstatus, 
         new_vital_status, os_event, 
         chemotherapy_drug = old_chemotherapy_drug_1, 
         hormone_therapy_drug = old_hormone_therapy_drug_1,
         days_from_treatment_start_to_sample_collection,
         days_from_drugs_end_to_sample_collection,
         chemo_length_indays, radiation_length_indays, hormone_length_indays,
         primary_site_desc,
         histology_desc = histology_desc.x, laterality_desc, 
         tnm_stage = tnm_stage1,
         ER_PR_status, ER_PR_HER_status,
         received_gcsf, gcsf_type,
         neutropenia_at_anytime)

write_csv(data_for_sharing, "Sequenced breast patients data 06272023.csv")
  


# End create variables
