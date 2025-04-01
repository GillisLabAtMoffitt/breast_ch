# Import library
library(tidyverse)
library(lubridate)

############################################################ I ### Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read.csv(paste0(#path_save,
  here::here(),
  "/processed data",
  "/Identified cleaned breast data_2024-12-17.csv")) %>% 
  mutate(mrn = as.character(mrn))
sequenced_patient_data <- read_rds(paste0(
  here::here(),
  "/processed data",
  "/Identified cleaned breast data_2024-12-17.rds"))

chart_dat <- 
  read_csv(paste0(#path_raw, 
    here::here(),
    "/chart_reviewed/cleaned chart review data.csv"))
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
chart_rev_cbc <- 
  readxl::read_xlsx(paste0(path_raw, "/chart_reviewed",
                           "/missing_cbc_2025-01-24.xlsx"),
                    na = c("NA")) %>% 
  janitor::clean_names()
Demographic <- 
  readxl::read_xlsx(paste0(path_raw, "/raw_data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "Demographics", na = c("Missing", "Unknown")) %>% 
  janitor::clean_names()
gcsf <- 
  readxl::read_xlsx(paste0(path_raw, "/raw_data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "Medication") %>% 
  janitor::clean_names()


############################################################ II ### G-CSF----
gcsf <- gcsf %>% 
  select(patient_id, drug_catalog_nm, drug_order_start_dtm) %>% 
  distinct() %>% 
  group_by(patient_id, drug_order_start_dtm) %>%
  summarise_at(vars(drug_catalog_nm), 
               str_c, collapse = "+") %>% 
  ungroup() %>% 
  # add mrn
  left_join(., Demographic %>% 
              select(mrn, patient_id), by = c("patient_id")) %>% 
  select(mrn, patient_id, 
         gcsf_type = drug_catalog_nm,
         gcsf_date = drug_order_start_dtm) %>% 
  # filter date before treatment
  inner_join(., sequenced_patient_data %>% 
               distinct(mrn, chemotherapy_start_date_1, chart_reviewed_chemotherapy_end_date_1, .keep_all = FALSE), 
             by = "mrn") %>% 
  filter(gcsf_date > chemotherapy_start_date_1 &
           gcsf_date < chart_reviewed_chemotherapy_end_date_1) %>% 
  # summarize
  group_by(mrn) %>%
  summarise_at(vars(gcsf_type, gcsf_date), 
               str_c, collapse = ", next regimen=") %>% 
  ungroup() %>% 
  mutate(received_gcsf = "Received G-CSF", .after = "mrn")
  
sequenced_patient_data <- sequenced_patient_data %>% 
  # Add G-CSF
  left_join(., gcsf, by = "mrn")
  
############################################################ III ### At pre-sample CBC----
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

# Clean baseline (pre-sample) CBC abstracted by Nancy
chart_dat <- chart_dat %>% 
  mutate(mrn = as.character(mrn)) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, date_of_cbc, wbc_k_u_l,
         rbc_mil_u_l, hemoglobin_g_d_l,
         mcv_fl, rdw_fl, platelets_k_u_l, 
         neutophils_k_u_l)

chart_rev_cbc <- chart_rev_cbc %>% 
  mutate(mrn = as.character(mrn)) %>% 
  filter(missing_datatype != "rdw") %>% 
  `colnames<-`(c("mrn", paste0(colnames(.)[2:ncol(.)], "_2nd_chart_rev"))) %>% 
  # mutate(specimen_collection_date_2nd_chart_rev = as.Date(specimen_collection_date_2nd_chart_rev))
  # mutate(across(contains("date"), ~ as.Date(.)))
  # mutate(across(contains("date"), ~ as.Date(as.numeric(.),
  #                                           origin = "1899-12-30")))
  mutate(across(contains("date"), ~ case_when(
    str_detect(., "T")                  ~ as.Date(., format = c("%Y-%m-%d")),
    TRUE                                ~ as.Date(as.numeric(.),
                                                 origin = "1899-12-30")
  )))


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
           lab_nm == "MCV" | lab_nm == "RDW" | lab_nm == "RBC") %>% 
  mutate(lab_nm = case_when(
    lab_nm == "Hemoglobin"            ~ "hemoglobin",
    lab_nm == "Platelet Count(k/uL)"  ~ "platelet",
    lab_nm == "WBC(k/uL)"             ~ "wbc",
    lab_nm == "MCV"                   ~ "mcv",
    lab_nm == "RDW"                   ~ "rdw",
    lab_nm == "RBC"                   ~ "rbc"
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
  # select neutrophil poly and bands
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
  select(mrn, lab_result_neutrophil = lab_result, lab_unit_neutrophil = lab_unit, lab_date_neutrophil = order_dtm) %>% 
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
  mutate(overall_neutrophil_lab_result = coalesce(lab_result_neutrophil, lab_result_neutrophlil_poly_band)) %>% 
  mutate(overall_neutrophil_lab_date = coalesce(lab_date_neutrophil, neutrophlil_poly_band_lab_date)) %>% 
  # I verify all units of selected values are k/uL
  mutate(overall_neutrophil_lab_units = "k/uL")

sequenced_patient_data <- 
  # Join CBC at pre-sample----
  full_join(cbc1, neutrophil, by = "mrn") %>% 
  full_join(., chart_dat, by = "mrn") %>% 
  full_join(., chart_rev_cbc, by = "mrn") %>% 
  # Add lab results
  mutate(lab_result_hemoglobin = coalesce(lab_result_hemoglobin, hemoglobin_g_d_l, lab_result_hemoglobin_at_presample_2nd_chart_rev)) %>% 
  mutate(lab_result_mcv = coalesce(lab_result_mcv, mcv_fl, lab_result_mcv_at_presample_2nd_chart_rev)) %>% 
  mutate(lab_result_platelet = coalesce(lab_result_platelet, platelets_k_u_l, lab_result_platelet_at_presample_2nd_chart_rev)) %>% 
  mutate(lab_result_rbc = coalesce(lab_result_rbc, rbc_mil_u_l, lab_result_rbc_at_presample_2nd_chart_rev)) %>% 
  mutate(lab_result_rdw = coalesce(lab_result_rdw, rdw_fl, lab_result_rdw_at_presample_2nd_chart_rev)) %>% 
  mutate(lab_result_wbc = coalesce(lab_result_wbc, wbc_k_u_l, lab_result_wbc_at_presample_2nd_chart_rev)) %>% 
  mutate(overall_neutrophil_lab_result = coalesce(overall_neutrophil_lab_result, neutophils_k_u_l, lab_result_neutrophil_at_presample_2nd_chart_rev)) %>% 
  # Add lab date
  mutate(lab_date_hemoglobin = case_when(
    !is.na(lab_result_hemoglobin)                   ~ coalesce(lab_date_hemoglobin, date_of_cbc, lab_date_hemoglobin_at_presample_2nd_chart_rev))) %>% 
  mutate(lab_date_mcv = case_when(
    !is.na(lab_result_mcv)                          ~ coalesce(lab_date_mcv, date_of_cbc, lab_date_mcv_at_presample_2nd_chart_rev))) %>% 
  mutate(lab_date_platelet = case_when(
    !is.na(lab_result_platelet)                     ~ coalesce(lab_date_platelet, date_of_cbc, lab_date_platelet_at_presample_2nd_chart_rev))) %>% 
  mutate(lab_date_rbc = case_when(
    !is.na(lab_result_rbc)                          ~ coalesce(lab_date_rbc, date_of_cbc, lab_date_rbc_at_presample_2nd_chart_rev))) %>% 
  mutate(lab_date_rdw = case_when(
    !is.na(lab_result_rdw)                          ~ coalesce(lab_date_rdw, date_of_cbc, lab_date_rdw_at_presample_2nd_chart_rev))) %>% 
  mutate(lab_date_wbc = case_when(
    !is.na(lab_result_wbc)                          ~ coalesce(lab_date_wbc, date_of_cbc, lab_date_wbc_at_presample_2nd_chart_rev))) %>% 
  mutate(overall_neutrophil_lab_date = case_when(
    !is.na(overall_neutrophil_lab_result)           ~ coalesce(overall_neutrophil_lab_date, date_of_cbc, lab_date_neutrophil_at_presample_2nd_chart_rev))) %>% 
  # Add lab unit
  mutate(lab_unit_hemoglobin = case_when(
    !is.na(lab_result_hemoglobin)                   ~ coalesce(lab_unit_hemoglobin, "g/dL"))) %>% 
  mutate(lab_unit_mcv = case_when(
    !is.na(lab_result_mcv)                          ~ coalesce(lab_unit_mcv, "FL"))) %>% 
  mutate(lab_unit_platelet = case_when(
    !is.na(lab_result_platelet)                     ~ coalesce(lab_unit_platelet, "k/uL"))) %>% 
  mutate(lab_unit_rbc = case_when(
    !is.na(lab_result_rbc)                          ~ coalesce(lab_unit_rbc, "mil/uL"))) %>% 
  mutate(lab_unit_rdw = case_when(
    !is.na(lab_result_rdw)                          ~ coalesce(lab_unit_rdw, "FL"))) %>% 
  mutate(lab_unit_wbc = case_when(
    !is.na(lab_result_wbc)                          ~ coalesce(lab_unit_wbc, "k/uL"))) %>% 
  mutate(overall_neutrophil_lab_units = case_when(
    !is.na(overall_neutrophil_lab_result)           ~ coalesce(overall_neutrophil_lab_units, "g/dL"))) %>% 
  # Clean
  rename(date_of_cbc_chart_reviewed = date_of_cbc) %>% 
  select(-c(wbc_k_u_l : neutophils_k_u_l), -c(ends_with("_2nd_chart_rev"))) %>% 
  `colnames<-`(c("mrn", paste0(colnames(.)[2:ncol(.)], "_at_presample"))) %>% 
  # Add them to data
  left_join(sequenced_patient_data, ., by = "mrn") %>% 
  # Reorganize variables---
  select(mrn : neutropenia_at_anytime,
         received_gcsf, gcsf_type, gcsf_date,
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
  filter(is.na(overall_neutrophil_lab_result_at_presample) |
           is.na(lab_result_hemoglobin_at_presample) |
           is.na(lab_result_mcv_at_presample) |
           is.na(lab_result_platelet_at_presample) |
           is.na(lab_result_rdw_at_presample) |
           is.na(lab_result_wbc_at_presample)
  ) %>% 
  filter(sample_treatment_sequence == "pre") %>% 
  mutate(lab_date_neutrophil_band_at_presample = neutrophlil_poly_band_lab_date_at_presample) %>% 
  select(mrn,
         specimen_collection_date,
         date_of_first_corresponding_treatment,
         lab_result_hemoglobin_at_presample, lab_unit_hemoglobin_at_presample, lab_date_hemoglobin_at_presample, 
         lab_result_mcv_at_presample, lab_unit_mcv_at_presample, lab_date_mcv_at_presample, 
         lab_result_platelet_at_presample, lab_unit_platelet_at_presample, lab_date_platelet_at_presample, 
         lab_result_rbc_at_presample, lab_unit_rbc_at_presample, lab_date_rbc_at_presample, 
         lab_result_rdw_at_presample, lab_unit_rdw_at_presample, lab_date_rdw_at_presample, 
         lab_result_wbc_at_presample, lab_unit_wbc_at_presample, lab_date_wbc_at_presample, 
         lab_result_neutrophil_at_presample, lab_unit_neutrophil_at_presample, lab_date_neutrophil_at_presample, 
         lab_result_neutrophlil_poly_band_at_presample, neutrophlil_poly_band_lab_date_at_presample, 
         lab_result_neutrophil_poly_at_presample, lab_unit_neutrophil_poly_at_presample, lab_date_neutrophil_poly_at_presample, 
         lab_result_neutrophil_bands_at_presample, lab_unit_neutrophil_bands_at_presample, 
         lab_date_neutrophil_band_at_presample,
         overall_neutrophil_lab_result_at_presample, overall_neutrophil_lab_units_at_presample, overall_neutrophil_lab_date_at_presample
         # contains("lab_date"),
         # # date_of_cbc_chart_reviewed_at_presample,
         # ends_with("_at_presample")
         ) %>% 
  mutate(missing_data_neutro = case_when(is.na(overall_neutrophil_lab_result_at_presample) ~ "neutro")) %>% 
  mutate(missing_data_hgb = case_when(is.na(lab_result_hemoglobin_at_presample) ~ "hgb")) %>% 
  mutate(missing_data_mcv = case_when(is.na(lab_result_mcv_at_presample) ~ "mcv")) %>% 
  mutate(missing_data_plt = case_when(is.na(lab_result_platelet_at_presample) ~ "plt")) %>% 
  mutate(missing_data_rdw = case_when(is.na(lab_result_rdw_at_presample) ~ "rdw")) %>% 
  mutate(missing_data_wbc = case_when(is.na(lab_result_wbc_at_presample) ~ "wbc")) %>% 
  unite(missing_datatype, starts_with("missing_data_"), sep = ", ", na.rm = TRUE) %>% 
  select(mrn, missing_datatype, everything(), -c(starts_with("missing_data_"))) %>% 
  arrange(missing_datatype)



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


# END Create new variable - (with CBC)

