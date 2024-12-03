# Import library
library(tidyverse)
library(lubridate)

############################################################ I ### Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read.csv(paste0(#path_save,
  here::here(),
  "/processed data",
  "/Identified cleaned breast data_2024-12-03.csv")) %>% 
  mutate(mrn = as.character(mrn))
sequenced_patient_data <- read_rds(paste0(
  here::here(),
  "/processed data",
  "/Identified cleaned breast data_2024-12-03.rds"))

# wbc <- read.csv(paste0(#path_save,
#   here::here(),
#   "/processed data/Cleaned WBC data.csv"))
# hgb_plt <- read.csv(paste0(#path_save,
#   here::here(),
#   "/processed data/Cleaned HgB-Plt data.csv"))
neutrophil <- 
  read_csv(paste0(path_save, "/processed data/Neutrophil lab data.csv"))

path_raw <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "sequential_samples_hossein")
cbc <- 
  readxl::read_xlsx(paste0(path_raw, "/raw_data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "CBC") %>% 
  janitor::clean_names()
Demographic <- 
  readxl::read_xlsx(paste0(path_raw, "/raw_data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "Demographics", na = c("Missing", "Unknown")) %>% 
  janitor::clean_names()

############################################################ II ### At pre-sample CBC----
# Create a temp data with pre and post dates variable
# dat <- sequenced_patient_data %>% 
#   mutate(mrn = as.character(mrn)) %>% 
#   mutate(presample_date = case_when(
#     sample_treatment_sequence == "pre"                       ~ specimen_collection_date
#   )) %>% 
#   mutate(postsample_date = case_when(
#     sample_treatment_sequence == "post"                      ~ specimen_collection_date
#   )) %>% 
#   group_by(mrn) %>% 
#   fill(presample_date, postsample_date, .direction = "updown") %>% 
#   ungroup() %>% 
#   select(mrn, presample_date, postsample_date) %>% 
#   distinct(mrn, .keep_all = TRUE)

# CBC cleaning
cbc <- cbc %>% 
  select(patient_id, lab_nm, lab_result, lab_unit, lab_date) %>% 
  mutate(lab_result = as.numeric(lab_result)) %>% 
  #add mrn
  left_join(., Demographic %>% 
              select(mrn, patient_id), by = c("patient_id"))


# CBC at pre-sample must be before or within 7 days from treatment start date 
# (aka date_of_first_corresponding_treatment)
cbc1 <- cbc %>% 
  filter((lab_nm == "Hemoglobin" & lab_unit == "g/dL") | 
           (lab_nm == "Platelet Count(k/uL)" & lab_unit == "k/uL") |
           (lab_nm == "WBC(k/uL)" & lab_unit == "k/uL") |
           lab_nm == "MCV" | lab_nm == "RDW") %>% 
  mutate(lab_nm = case_when(
    lab_nm == "Hemoglobin"            ~ "hemoglobin",
    lab_nm == "Platelet Count(k/uL)"  ~ "platelet",
    lab_nm == "WBC(k/uL)"             ~ "wbc",
    lab_nm == "MCV"                   ~ "mcv",
    lab_nm == "RDW"                   ~ "rdw",
  )) %>% 
  # add date of pre-sample and select closest cbc
  inner_join(., sequenced_patient_data %>% 
               filter(sample_treatment_sequence == "pre") %>%
               distinct(mrn, date_of_first_corresponding_treatment, .keep_all = FALSE), 
             by = "mrn") %>% 
  # filter date before treatment
  filter(lab_date <= date_of_first_corresponding_treatment + days(7)) %>% 
  # pick the closest date to treatment
  mutate(int = abs(interval(start = lab_date, end = date_of_first_corresponding_treatment)/
                     duration(n = 1, units = "days"))) %>% 
  arrange(mrn, lab_nm, int) %>% 
  distinct(mrn, lab_nm, .keep_all = TRUE) %>% 
  # pivot
  select(mrn, lab_nm, lab_result, lab_unit, lab_date) %>% 
  pivot_wider(id_cols = mrn, 
              names_from = lab_nm, 
              values_from = c(lab_result, lab_unit, lab_date), 
              names_vary = "slowest")

# Neutro----
# I take neutro or auto first and when missing add the others...
neutrophil_other <- neutrophil %>% 
  mutate(mrn = as.character(mrn)) %>% 
  # select neutroplil poly and bands
  filter(lab_nm != "Neutrophil") %>% 
  select(mrn, lab_nm, lab_result, lab_unit, lab_date = order_dtm) %>% 
  mutate(lab_result = case_when(
    lab_unit == "k/uL"           ~ lab_result,
    lab_unit == "cells/uL"       ~ lab_result / 1000
  )) %>% 
  # add date of pre-sample and select closest cbc
  inner_join(., sequenced_patient_data %>% 
               filter(sample_treatment_sequence == "pre") %>%
               distinct(mrn, date_of_first_corresponding_treatment, .keep_all = FALSE), 
             by = "mrn") %>% 
  # filter date before treatment
  filter(lab_date <= date_of_first_corresponding_treatment + days(7)) %>% 
  # pick the closest date to treatment
  mutate(int = abs(interval(start = lab_date, end = date_of_first_corresponding_treatment)/
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
  rename(neutrophlil_poly_band_lab_date = lab_date_neutrophil_bands) %>% 
  filter(!is.na(lab_result_neutrophlil_poly_band))

neutrophil <- neutrophil %>% 
  mutate(mrn = as.character(mrn)) %>% 
  filter(lab_nm == "Neutrophil") %>% 
  select(mrn, lab_result_neutroplil = lab_result, lab_unit_neutroplil = lab_unit, lab_date_neutrophil = order_dtm) %>% 
  # add date of pre-sample and select closest cbc
  inner_join(., sequenced_patient_data %>% 
               filter(sample_treatment_sequence == "pre") %>%
               distinct(mrn, date_of_first_corresponding_treatment, .keep_all = FALSE), 
             by = "mrn") %>% 
  # filter date before treatment
  filter(lab_date_neutrophil <= date_of_first_corresponding_treatment + days(7)) %>% 
  # pick the closest date to treatment
  mutate(int = abs(interval(start = lab_date_neutrophil, end = date_of_first_corresponding_treatment)/
                     duration(n = 1, units = "days"))) %>% 
  arrange(mrn, int) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(-c(date_of_first_corresponding_treatment, int))

neutrophil <- full_join(neutrophil, neutrophil_other, by = "mrn") %>% 
  mutate(overall_neutroplil_lab_result = coalesce(lab_result_neutroplil, lab_result_neutrophlil_poly_band)) %>% 
  mutate(overall_neutroplil_lab_date = coalesce(lab_date_neutrophil, neutrophlil_poly_band_lab_date)) %>% 
  # I verify all units of selected values are k/uL
  mutate(overall_neutroplil_lab_units = "k/uL")

sequenced_patient_data <- 
  # Join CBC at pre-sample----
  full_join(cbc1, neutrophil, by = "mrn") %>% 
  `colnames<-`(c("mrn", paste0(colnames(.)[2:ncol(.)], "_at_presample"))) %>% 
  # Add them to data
  left_join(sequenced_patient_data, ., by = "mrn") %>% 
  # Reorganize variables---
  select(mrn : neutropenia_at_anytime,
         ends_with("_at_presample"),
         everything(), -neutropenia_at_anytime)


# Save clean data with new variables (with CBC variables)
write_csv(sequenced_patient_data, 
          paste0(here::here(), 
                 "/processed data",
                 "/Identified cleaned breast data with CBC at presample_",
                 today(), ".csv"))
write_rds(sequenced_patient_data, 
          paste0(here::here(), 
                 "/processed data",
                 "/Identified cleaned breast data with CBC at presample_",
                 today(), ".rds"))
write_csv(sequenced_patient_data, 
          paste0(path_save, 
                 "/processed data",
                 "/Identified cleaned breast data with CBC at presample_",
                 today(), ".csv"))

missing_cbc <- sequenced_patient_data %>% 
  filter(is.na(overall_neutroplil_lab_result_at_presample) |
           is.na(lab_result_hemoglobin_at_presample) |
           is.na(lab_result_mcv_at_presample) |
           is.na(lab_result_platelet_at_presample) |
           is.na(lab_result_rdw_at_presample) |
           is.na(lab_result_wbc_at_presample)
  ) %>% 
  filter(sample_treatment_sequence == "pre") %>% 
  select(mrn,
         ends_with("_at_presample"))
write_csv(missing_cbc, 
          paste0(here::here(), 
                 "/processed data",
                 "/missing_cbc_",
                 today(), ".csv"))
write_csv(missing_cbc, 
          paste0(path_save, 
                 "/processed data",
                 "/missing_cbc_",
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
                 "/De-identified cleaned breast data with CBC at presample_", 
                 today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data", 
                 "/De-identified cleaned breast data with CBC at presample_", 
                 today(), ".csv"))


# Cytopenia----

























# END Create new variable - (with CBC)

