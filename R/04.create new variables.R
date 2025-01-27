# Import library
library(tidyverse)
library(lubridate)

############################################################ I ### Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read_rds(paste0(
  here::here(),
  "/processed data",
  "/Identified breast data with re-classified sequenced sequential sample and chart review primary clinical_2024-12-17.rds"))

subsequent_cancer <-
  read_csv(paste0(here::here(), "/processed data/patient subsequent cancer after breast cancer - hossein.csv"))


############################################################ II ### Adding de-identified and extra variables----
# summarize subsequent cancer
subsequent_cancer <- subsequent_cancer %>% 
  mutate(mrn = as.character(mrn)) %>% 
  rename(dx_date = dx_dt) %>% 
  group_by(mrn) %>% 
  summarise_at(vars(dx_date, histology_cd, histology_desc,
                    primary_site_cd, primary_site_group_desc), 
               str_c, collapse = "; ") %>%
  ungroup() %>% 
  `colnames<-`(c("mrn", paste0("subsequent_cancer_", colnames(.)[2:ncol(.)])))

sequenced_patient_data <- sequenced_patient_data %>%
  # add subsequent cancer----
  left_join(., subsequent_cancer, by = "mrn") %>% 
  # add new variables----
  # create age
  mutate(age_at_diagnosis = round(interval(start = date_of_birth, end = date_of_diagnosis)/
                                    duration(n = 1, units = "years"), 1)
  ) %>%
  mutate(age_at_sample = round(interval(start = date_of_birth, end = specimen_collection_date)/
                                 duration(n = 1, units = "years"), 1)
  ) %>% 
  mutate(year_at_sample = year(specimen_collection_date)) %>% 
  # treatment var
  mutate(date_of_first_corresponding_treatment = case_when(
    treatment_type == "chemo"              ~ chemotherapy_start_date_1,
    treatment_type == "chemorad"           ~ chemotherapy_start_date_1,
    treatment_type == "hormone"            ~ hormone_therapy_start_date_1,
    treatment_type == "radiation"          ~ radiation_start_date_1
  )) %>% 
  mutate(age_at_first_treatment = round(interval(start = date_of_birth, end = date_of_first_corresponding_treatment)/
                                          duration(n = 1, units = "years"), 1)
  ) %>% 
  mutate(age_at_first_chemo = round(interval(start = date_of_birth, end = chemotherapy_start_date_1)/
                                      duration(n = 1, units = "years"), 1)
  ) %>% 
  mutate(age_at_first_hormone = round(interval(start = date_of_birth, end = hormone_therapy_start_date_1)/
                                        duration(n = 1, units = "years"), 1)
  ) %>% 
  mutate(age_at_first_radiation = round(interval(start = date_of_birth, end = radiation_start_date_1)/
                                          duration(n = 1, units = "years"), 1)
  ) %>% 
  mutate(days_from_treatment_start_to_sample_collection = case_when(
    treatment_type == "chemo"              ~ interval(start = chemotherapy_start_date_1,
                                                      end = specimen_collection_date)/
      duration(n = 1, units = "days"),
    treatment_type == "chemorad"           ~ interval(start = chemotherapy_start_date_1,
                                                      end = specimen_collection_date)/
      duration(n = 1, units = "days"),
    treatment_type == "hormone"            ~ interval(start = hormone_therapy_start_date_1,
                                                      end = specimen_collection_date)/
      duration(n = 1, units = "days"),
    treatment_type == "radiation"          ~ interval(start = radiation_start_date_1,
                                                      end = specimen_collection_date)/
      duration(n = 1, units = "days")
  )) %>% 
  mutate(days_from_chemo_start_to_radiation = case_when(
    treatment_type == "chemorad"           ~ interval(start = chemotherapy_start_date_1,
                                                      end = radiation_start_date_1)/
      duration(n = 1, units = "days")
  )) %>%
  mutate(chemo_length_indays = case_when(
    treatment_type == "chemo" |
      treatment_type == "chemorad"         ~ interval(start = chemotherapy_start_date_1,
                                                      end = chart_reviewed_chemotherapy_end_date_1)/
      duration(n = 1, units = "days")
  )) %>% 
  mutate(radiation_length_indays = case_when(
    treatment_type == "chemorad" |
      treatment_type == "radiation"        ~ interval(start = radiation_start_date_1,
                                                      end = radiation_end_date1_1)/
      duration(n = 1, units = "days")
  )) %>% 
  mutate(hormone_length_indays = case_when(
    treatment_type == "hormone"            ~ interval(start = hormone_therapy_start_date_1,
                                                      end = hormone_therapy_end_date1_1)/
      duration(n = 1, units = "days")
  )) %>% 
  # Survival----
  # OS
  mutate(os_event = case_when(
    dead_or_alive == "Alive"            ~ 0,
    dead_or_alive == "Dead"             ~ 1
  ), .after = dead_or_alive) %>% 
  mutate(os_time_from_dx_months = interval(start = date_of_diagnosis,
                                           end = date_of_last_followup)/
           duration(n = 1, unit = "months"),
         .after = os_event) %>%
  mutate(os_time_from_tx_start_months = interval(start = date_of_first_corresponding_treatment,
                                                 end = date_of_last_followup)/
           duration(n = 1, unit = "months"),
         .after = os_event) %>%
  # PFS
  mutate(pfs_event = case_when(
    is.na(progression_type)             ~ 0,
    !is.na(progression_type)            ~ 1
  ), .after = progression_type) %>% 
  mutate(pfs_date = coalesce(date_of_progression, date_of_last_followup), 
         .after = pfs_event) %>% 
  mutate(pfs_time_months = interval(start = date_of_first_corresponding_treatment, end = pfs_date)/
           duration(n = 1, unit = "months"),
         .after = pfs_event)
  
# Reorganize variables---
sequenced_patient_data <- sequenced_patient_data %>%
  select(mrn, deidentified_patient_id, party_id,
         # Sample
         sample_id, sample_family_id, 
         submitted_ID_for_added_samples, received_ID_for_added_samples,
         Sample_Name_Fastq_file, sample_name_secondary,
         specimen_collection_date, age_at_sample, year_at_sample,
         # Treatment
         treatment_type, discrepancies_old_new_sample_sequence,
         sample_treatment_sequence = time_to_treatment, 
         sequence_sample_vs_treatment, treatment_received,
         date_of_first_corresponding_treatment,
         chart_reviewed_end_date_of_first_corresponding_treatment,
         age_at_first_treatment, age_at_first_chemo,
         age_at_first_hormone, age_at_first_radiation,
         days_from_treatment_start_to_sample_collection,
         
         # interval_presample_chemo : blood_bf_chemo_rad, had_good_sample_chemo_rad,
         
         chemotherapy_start_date_1 : first_chemo_rad_date,
         chemo_length_indays, radiation_length_indays,
         days_from_chemo_start_to_radiation,
         hormone_length_indays, 
         chemotherapy_start_date_1, chemotherapy_drug_1, chart_reviewed_chemotherapy_end_date_1,
         hormone_therapy_start_date_1, hormone_therapy_drug_1, hormone_therapy_end_date1_1,
         radiation_start_date_1, type_of_radiation_1, total_dose_radiation_1, 
         boost_dose_c_gy_radiation_1, radiation_end_date1_1, 
         first_chemo_rad_date,
         chemotherapy_drug_2 : radiation_end_date1_5,
         # Cancer
         age_at_diagnosis, 
         date_of_diagnosis : tnm_cs_mixed_group_stage4,
         primary_site_cd, primary_site_desc, 
         histology_code, histology, laterality,
         tnm_stage, clinical_tnm_group_stage, 
         tnm_cs_mixed_group_stage = tnm_cs_mixed_group_stage1,
         class_of_case_cd,
         starts_with("subsequent"),
         er : oncotype_dx,
         ln_positive_axillary_level_i_ii,
         contains("multigene_signature"),
         contains("ssf"),
         # received_gcsf, gcsf_type, 
         neutropenia_at_anytime,
         # Demo
         date_of_birth,
         gender_cancer_registry : gender_cerner,
         gender, race, ethnicity, 
         csv_ethnicity_cancer_registry : race_cerner,
         smoking_status, smoking_amount : quit_date,
         marital_status_current,
         starts_with("csv_what_is_the_last_grade"),
         starts_with("csv_what_is_your_current_occupation"),
         starts_with("csv_what_was_your_occupation"),
         # Survival
         os_time_from_dx_months,
         os_time_from_tx_start_months,
         os_event, dead_or_alive, date_of_last_followup,
         metastatic_site_at_diagnosis_desc, 
         pfs_event : date_of_progression, progression_type,
         contains("response_to_neoadjuvant_therapy"),
         # type_of_first_recurrence_desc,
         starts_with("csv_cause_of_death"),
         max_date_of_last_contact_or_death_cancer_registry :
           csv_vital_status_cerner,
         
         everything(), -new_last_date
  )

# Save clean data with new variables (without CBC variables)
write_csv(sequenced_patient_data, 
          paste0(here::here(), 
                 "/processed data",
                 "/Identified cleaned breast data_",
                 today(), ".csv"))
write_rds(sequenced_patient_data, 
          paste0(here::here(), 
                 "/processed data",
                 "/Identified cleaned breast data_",
                 today(), ".rds"))
write_csv(sequenced_patient_data, 
          paste0(path_save, 
                 "/processed data",
                 "/Identified cleaned breast data_",
                 today(), ".csv"))


# Save de-identified data
deids_data <- sequenced_patient_data %>% 
  select(-c(mrn, sample_id, sample_family_id, sample_name_secondary, party_id,
            submitted_ID_for_added_samples, received_ID_for_added_samples,
            Sample_Name_Fastq_file,
            contains("date"),
            blood_bf_chemo_rad : sequenced_samples
  ))

write_csv(deids_data, 
          paste0(path_save, 
                 "/processed data", 
                 "/De-identified cleaned breast data_", 
                 today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data", 
                 "/De-identified cleaned breast data_", 
                 today(), ".csv"))


# END Create new variable - (without CBC)
