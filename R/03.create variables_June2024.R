# Import library
library(tidyverse)
library(lubridate)

############################################################ I ### Load data
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

blood_patients <- 
  read_rds(paste0(here::here(), "/blood_patients2.rds"))

# blood_patients3 <-
#   read_rds(paste0(here::here(), "/blood_patients2_3years.rds"))

later_blood_patients <-
  read_rds(paste0(here::here(), "/blood_patients_06272023.rds"))

# later_blood_patients2 <-
#   read_rds(paste0(here::here(), "/blood_patients02282023.rds"))


CH_status <- 
  readxl::read_xlsx(paste0(
    here::here(), "/list_to_send_out_BickWHO_05.14.2024.xlsx"), na = "NA")
# The core swapped some samples with other samples from the same date 
# that are not in the raw file
# Need the submission form and received ids data
received_samples <- 
  readxl::read_xlsx(paste0(
    here::here(), "/Submission Form-CICPT-3757- Gillis_AgilentDNA_Panel_08.25.23.xlsx"), 
    sheet = "Submission Form", skip = 9)
submitted_samples <- 
  readxl::read_xlsx(paste0(
    here::here(), "/Submission Form-CICPT-3757- Gillis_AgilentDNA_Panel_08.25.23.xlsx"), 
    sheet = "Sample List")

# Same for the first sequencing data (only 1 sample)
fixed_ids_1 <- 
  read_csv(paste0(here::here(), "/breast sample swapped ids in the first sequencing run_CCL.csv"))
# samples <- read_rds(paste0(here::here(), "/samples sequenced as of june 2023.rds"))
# 
# old_blood_patients <- 
#   read_rds("/Users/colinccm/Documents/GitHub/Gillis/breast_ch/blood_patients.rds")
# path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Breast CH")
# missing_data <- 
#   read_csv(paste0(path, "/raw data/Missing data.csv"))

subsequent_cancer <- 
  read_csv(paste0(here::here(), "/patient subsequent cancer after breast cancer - hossein.csv"))

path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Breast CH")
cbc <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "CBC") %>% 
  janitor::clean_names()
Demographic <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "Demographics", na = c("Missing", "Unknown")) %>% 
  janitor::clean_names()
neutrophlil <- 
  read_csv(paste0(path, "/processed data/Neutrophil lab data.csv"))


############################################################ II ### Merge data and limit to patients sequenced
# pivot longer sample id
blood_patients <- blood_patients %>% 
  # Need 1 id per row to merge with sample sequenced
  mutate(sample_count = sapply(str_split(sample_id, ";"), length)) %>% 
  separate_wider_delim(cols = sample_id, delim = "; ",
                       names = c(paste("sample_id", 1:max(.$sample_count), sep = "_")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  pivot_longer(cols = c(starts_with("sample_id_")), 
               names_to = NULL, values_to = "sample_id", 
               values_drop_na = TRUE) %>% 
  select(mrn, 
         sample_family_id, sample_id,
         everything(), -sample_count) 

# Fix the swapped sample ids 
# I verified that the collection dates are the same
fixed_ids <- submitted_samples %>% 
  rename(submitted_ID = SampleID) %>% 
  full_join(., received_samples %>% 
              filter(!is.na(`LV Alias DNA ID`)) %>% 
              rename(received_ID = `DNA-Aliquot-LV`), 
            by = c("Sample_Name" = "Fastq file-name")) %>% 
  filter(!is.na(received_ID)) %>% 
  bind_rows(., fixed_ids_1) %>% 
  mutate(MRN = as.character(MRN)) %>% 
  select(MRN, deidentified_patient_id, 
         submitted_ID, received_ID, Sample_Name)

CH_status <- CH_status %>% 
  # Cleanup
  select(-c(patient_id, ids)) %>%
  mutate(sequenced_samples = "Yes") %>% 
  # add the swapped sample ids
  full_join(., fixed_ids %>% 
              mutate(`Sample Name (Secondary)` = received_ID), 
            by = c("MRN", "Sample Name (Secondary)", 
                   "Sample_Name_Fastq_file" = "Sample_Name")) %>% 
  select(MRN, Sample_Name_Fastq_file,
         "Sample Name (Secondary)", 
         submitted_ID_for_added_samples = submitted_ID, 
         received_ID_for_added_samples = received_ID,
         `Collection Date` : sequenced_samples) %>% 
  # Fill the treatment for the extra samples that Nancy choose later
  group_by(MRN) %>% 
  fill(treatment, .direction = "updown") %>% 
  ungroup() %>% 
  # mutate(CH_status = case_when(
  #   !is.na(VARIANT_C) &
  #     !is.na(VARIANT_P) &
  #     !is.na(VAF) &
  #     !is.na(DEPTH)             ~ "CH",
  #   TRUE                        ~ "No CH"
  # ))
  distinct(MRN, Sample_Name_Fastq_file, .keep_all = TRUE) %>% 
  select(-c(VARIANT_C : DEPTH))

sequenced_patients <- blood_patients %>% 
  # mutate(has_clinical = "Yes") %>% 
  full_join(., CH_status %>% 
              # mutate(has_CH = "Yes") %>% 
              mutate(merging_ids = 
                       coalesce(submitted_ID_for_added_samples, `Sample Name (Secondary)`)) %>% 
            mutate(merging_ids = str_to_lower(merging_ids)),
            by = c("mrn" = "MRN", "sample_id" = "merging_ids"
            )) %>% 
  mutate(deidentified_patient_id =
           str_match(Sample_Name_Fastq_file, "(.*_\\d+)_")[,2]) %>%
  # Limit to patients sequenced
  filter(!is.na(sequenced_samples)) %>%
  select(mrn, deidentified_patient_id, sample_id, sample_family_id,
         submitted_ID_for_added_samples, received_ID_for_added_samples,
         Sample_Name_Fastq_file,
         time_to_treatment = status, treatment_type = treatment, 
         interval_presample_chemo, time_chemo_seqsample, interval_prepost_sample_chemo,
         interval_presample_rad, time_rad_seqsample, interval_prepost_sample_rad,
         interval_presample_hormone, time_hormone_seqsample, interval_prepost_sample_hormone,
         # interval_presample_chemrad, time_chemorad_seqsample, interval_prepost_sample_chemorad, interval_chemo_rad, 
         interval_presample_chemo_rad, time_chemo_1st_seqsample, 
         time_1st_seqsample_rad, time_rad_2nd_seqsample, interval_prepost_sample_chemo_rad,
         everything(), sequenced_samples) %>% 
  # Change time pre-samples to treatment to be negative values as Nancy's request
  mutate(across(c(interval_presample_chemo, interval_presample_rad, 
                  interval_presample_hormone, interval_presample_chemrad, 
                  interval_presample_chemo_rad), ~ -.))


# Add the variables created in the last clinical version data of June 
# but doesn't have the selected samples info
sequenced_patients <- later_blood_patients %>% 
  # mutate(has_later_clin = "Yes") %>% 
  arrange(mrn, date_of_diagnosis1) %>% 
  select(mrn, birth_dt,
         gender, race, ethnicity, 
         primary_site_cd, primary_site_desc, 
         histology_code = histology, laterality = laterality_desc,
         tnm_stage = tnm_stage1,
         class_of_case_cd,
         metastatic_site_at_diagnosis_desc, 
         type_of_first_recurrence_desc,
         new_vital_status, new_last_date, smoking_status,
         received_gcsf, gcsf_type, neutropenia_at_anytime,
         ssf1_nm : progesterone_receptor_total_allred_score_desc,
         ER_results, PR_results, HER_results, 
         ER_PR_status, ER_PR_HER_status
         ) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  left_join(sequenced_patients, ., by = c("mrn", "date_of_birth" = "birth_dt")) %>% 
  mutate(race = coalesce(race, race_cancer_registry_1)) %>% 
  mutate(ethnicity = coalesce(ethnicity, ethnicity_cerner)) %>% 
  select(-age_at_diagnosis)

write_csv(sequenced_patients, "Sequenced breast data with sequential sample info and primary clinical.csv")


# Prep CBC for CHRS calculation----
# cbc1 <- cbc
cbc <- cbc1



hb_plt <- cbc %>% 
  filter((lab_nm == "Hemoglobin" & lab_unit == "g/dL") | 
           lab_nm == "Platelet Count(k/uL)" & lab_unit == "k/uL") %>% 
  mutate(lab_result = as.numeric(lab_result)) %>% 
  #add mrn
  left_join(., Demographic %>% 
              select(mrn, patient_id), by = c("patient_id")) 
write_csv(hb_plt, paste0(here::here(), "/processed data/Cleaned HgB-Plt data.csv"))
write_csv(hb_plt, paste0(path_save, "/processed data/Cleaned HgB-Plt data.csv"))


hb_plt <- hb_plt %>% 
  # add date of pre-sample and select closest cbc
  inner_join(., sequenced_patients %>% 
              filter(time_to_treatment == "pre") %>%
              distinct(mrn, specimen_collection_date, .keep_all = FALSE), ############# remove when change CH
             by = "mrn") %>% 
  mutate(int = abs(interval(start = lab_date, end = specimen_collection_date)/
           duration(n = 1, units = "days"))) %>% 
  arrange(mrn, lab_nm, int) %>% 
  distinct(mrn, lab_nm, .keep_all = TRUE) %>% 
  # pivot
  select(mrn, lab_nm, lab_result, lab_unit, lab_date) %>% 
  pivot_wider(id_cols = mrn, 
              names_from = lab_nm, 
              values_from = c(lab_result, lab_unit, lab_date), 
              names_vary = "slowest")

# I take neutro or auto first and when missing add the others...
neutrophlil_other <- neutrophlil %>% 
  mutate(mrn = as.character(mrn)) %>% 
  # select neutroplil poly and bands
  filter(lab_nm != "Neutrophil") %>% 
  select(mrn, lab_nm, lab_result, lab_unit, lab_date = order_dtm) %>% 
  mutate(lab_result = case_when(
    lab_unit == "k/uL"           ~ lab_result,
    lab_unit == "cells/uL"       ~ lab_result / 1000
  )) %>% 
  # add date of pre-sample and select closest cbc
  inner_join(., sequenced_patients %>% 
               filter(time_to_treatment == "pre") %>%
               select(mrn, specimen_collection_date), by = "mrn") %>%
  # select the closest value to pre-sample
  mutate(int = abs(interval(start = lab_date, end = specimen_collection_date)/
                     duration(n = 1, units = "days"))) %>% 
  arrange(mrn, lab_nm, int) %>% 
  distinct(mrn, lab_nm, .keep_all = TRUE) %>% 
  # pivot and add ploy + bands values together
  mutate(lab_nm = str_replace(lab_nm, " ", "_"),
         lab_nm = str_to_lower(lab_nm)) %>% 
  pivot_wider(id_cols = mrn, 
              names_from = lab_nm, 
              values_from = c(lab_result, lab_unit, lab_date), 
              names_vary = "slowest") %>% 
  mutate(lab_result_neutrophlil_poly_band = lab_result_neutrophil_bands + lab_result_neutrophil_poly) %>% 
  mutate(neutrophlil_poly_band_lab_date = lab_date_neutrophil_bands) %>% 
  filter(!is.na(lab_result_neutrophlil_poly_band))

neutrophlil <- neutrophlil %>% 
  mutate(mrn = as.character(mrn)) %>% 
  filter(lab_nm == "Neutrophil") %>% 
  select(mrn, lab_result_neutroplil = lab_result, lab_unit_neutroplil = lab_unit, lab_date_neutroplil = order_dtm) %>% 
  # add date of pre-sample and select closest cbc
  inner_join(., sequenced_patients %>% 
               filter(time_to_treatment == "pre") %>%
               select(mrn, specimen_collection_date), by = "mrn") %>% 
  # select the closest value to pre-sample
  mutate(int = abs(interval(start = lab_date_neutroplil, end = specimen_collection_date)/
                     duration(n = 1, units = "days"))) %>% 
  arrange(mrn, int) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(-c(specimen_collection_date, int))

neutrophlil <- full_join(neutrophlil, neutrophlil_other, by = "mrn") %>% 
  mutate(overall_neutroplil_lab_result = coalesce(lab_result_neutroplil, lab_result_neutrophlil_poly_band)) %>% 
  mutate(overall_neutroplil_lab_date = coalesce(lab_date_neutroplil, neutrophlil_poly_band_lab_date)) %>% 
  # I verify all units of selected values are k/uL
  mutate(overall_neutroplil_lab_units = "k/uL")
write_csv(neutrophlil, paste0(here::here(), "/processed data/Cleaned ANC data.csv"))
write_csv(neutrophlil, paste0(path_save, "/processed data/Cleaned ANC data.csv"))

# Cytopenias were defined by using World Health Organization criteria ######### Do it later
# anemia = hemoglobin concentration <13.0 g/dl in male participants 
# and <12.0 g/dl in female participants; 
# thrombocytopenia = platelet counts <150 × 109 cells/l; 
# and neutropenia = absolute neutrophil count <1.8 × 109 cells/l

# cells/ l == cells / 1000ul == cells /ul
# 
# 1000cells /ul == 1000*1000cells / 1000ul == 1000+1000cells / L == 1 x 106

# WBC at pre-sample----
wbc <- cbc %>% 
  # mutate(mrn = as.character(mrn)) %>% 
  filter(lab_nm == "WBC(k/uL)" & lab_unit == "k/uL") %>% 
  mutate(lab_result = as.numeric(lab_result)) %>% 
  #add mrn
  left_join(., Demographic %>% 
              select(mrn, patient_id), by = c("patient_id")) %>% 
  select(mrn, lab_result_wbc = lab_result, lab_unit_wbc = lab_unit, lab_date_wbc = order_dtm)
write_csv(wbc, paste0(here::here(), "/processed data/Cleaned WBC data.csv"))
write_csv(wbc, paste0(path_save, "/processed data/Cleaned WBC data.csv"))

wbc <- wbc %>% 
  # add date of pre-sample and select closest cbc
  inner_join(., sequenced_patients %>% 
               filter(time_to_treatment == "pre") %>%
               select(mrn, specimen_collection_date), by = "mrn") %>% 
  # select the closest value to pre-sample
  mutate(int = abs(interval(start = lab_date_wbc, end = specimen_collection_date)/
                     duration(n = 1, units = "days"))) %>% 
  arrange(mrn, int) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(-c(specimen_collection_date, int))

# CBC at pre-sample----
cbc_at_presample <- full_join(wbc, hb_plt, by = "mrn") %>% 
  full_join(., neutrophlil, by = "mrn")


# summarize subsequent cancer----
subsequent_cancer <- subsequent_cancer %>% 
  mutate(mrn = as.character(mrn)) %>% 
  group_by(mrn) %>% 
  summarise_at(vars(dx_dt, histology_cd, histology_desc,
                    primary_site_cd, primary_site_group_desc), 
               str_c, collapse = "; ") %>%
  ungroup() %>% 
  `colnames<-`(c("mrn", paste0("subsequent_cancer_", colnames(.)[2:ncol(.)])))

sequenced_patients <- sequenced_patients %>%
  # add subsequent cancer
  left_join(., subsequent_cancer, by = "mrn") %>% 
  # create age
  mutate(age_at_diagnosis = round(interval(start = date_of_birth, end = date_of_diagnosis1)/
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
    treatment_type == "chemorad"           ~ interval(start = chemotherapy_end_date1_1,
                                                             end = radiation_start_date_1)/
      duration(n = 1, units = "days")
  )) %>%
  mutate(chemo_length_indays = case_when(
    treatment_type == "chemo" |
      treatment_type == "chemorad"         ~ interval(start = chemotherapy_start_date_1,
                                                      end = chemotherapy_end_date1_1)/
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
  # Survival
  mutate(os_time_from_dx_months = interval(start = date_of_diagnosis1, 
                                           end = new_last_date)/
           duration(n = 1, units = "months")) %>% 
  mutate(os_time_from_tx_start_months = interval(start = date_of_first_corresponding_treatment,
                                                 end = new_last_date)/
           duration(n = 1, units = "months")) %>% 
  mutate(os_event = case_when(
    new_vital_status == "ALIVE"                 ~ 0,
    new_vital_status == "DEAD"                  ~ 1
  )) %>% 
  # CBC
  left_join(., cbc_at_presample, by = "mrn") %>% 
  select(mrn, deidentified_patient_id, sample_id, sample_family_id, 
         submitted_ID_for_added_samples, received_ID_for_added_samples,
         Sample_Name_Fastq_file,
         `Sample Name (Secondary)`,
         specimen_collection_date, age_at_sample, year_at_sample,
         sample_treatment_sequence = time_to_treatment, treatment_type, 
         date_of_first_corresponding_treatment,
         age_at_first_treatment, age_at_first_chemo,
         age_at_first_hormone, age_at_first_radiation,
         days_from_treatment_start_to_sample_collection,
         days_from_chemo_start_to_radiation,
         chemo_length_indays, radiation_length_indays,
         hormone_length_indays, 
         chemotherapy_start_date_1, chemotherapy_drug_1, chemotherapy_end_date1_1,
         hormone_therapy_start_date_1, hormone_therapy_drug_1, hormone_therapy_end_date1_1,
         radiation_start_date_1, boost_dose_c_gy_1, radiation_end_date1_1,
         
         os_time_from_dx_months,
         os_time_from_tx_start_months,
         os_event, new_vital_status,
         metastatic_site_at_diagnosis_desc, 
         type_of_first_recurrence_desc,
         starts_with("csv_cause_of_death"),
         
         gender, race, ethnicity, 
         date_of_birth, date_of_diagnosis = date_of_diagnosis1,
         age_at_diagnosis, 
         primary_site_cd, primary_site_desc, 
         histology_code, histology, laterality,
         tnm_stage, clinical_tnm_group_stage, tnm_cs_mixed_group_stage1,
         class_of_case_cd,
         starts_with("subsequent"),
         received_gcsf, gcsf_type, 
         
         lab_result_wbc : overall_neutroplil_lab_units,
         neutropenia_at_anytime,
         ER_results, PR_results, HER_results, 
         ER_PR_status, ER_PR_HER_status,
         ssf1_nm : progesterone_receptor_total_allred_score_desc,
         
         smoking_status, marital_status_current,
         starts_with("csv_what_is_the_last_grade"),
         starts_with("csv_what_is_your_current_occupation"),
         starts_with("csv_what_was_your_occupation"),
         
  )


write_csv(sequenced_patients, "Identified breast data with sequenced sequential sample and clinical_07092024.csv")

deids_data <- sequenced_patients %>% 
  select(-c(mrn, sample_id, sample_family_id, 
         submitted_ID_for_added_samples, received_ID_for_added_samples,
         Sample_Name_Fastq_file,
         contains("date")
))

write_csv(deids_data, "De-identified breast data for sequenced sequential sample_07092024.csv")
         