# Import library
library(tidyverse)
library(lubridate)

############################################################ I ### Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read_rds(paste0(
  here::here(),
  "/processed data",
  "/Identified cleaned breast data with CBC at presample_2025-01-27.rds"))

updated_end_dates <- read_rds(paste0(fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                                              "Breast_R01"),
                                     "/other_info",
                                     "/Sample_updated_date_with_interval_03.24.2025.rds"))

mrn_modif <- read_csv("mrn_hormone_end_date_is_na.csv")

chart_dat <- 
  read_csv(paste0(#path_raw, 
    here::here(),
    "/chart_reviewed/cleaned chart review data.csv"))

hgb_plt <- read.csv(paste0(#path_save,
  here::here(),
  "/processed data/Cleaned HgB-Plt data.csv"))

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


############################################################ II ### Prep new data date----
updated_end_dates <- updated_end_dates %>% 
  select(deidentified_patient_id, specimen_collection_date, 
         sample_treatment_sequence,
         sequence_sample_vs_treatment, 
         treatment_type, everything(), -day_between_pre_lastpost) %>% 
  `colnames<-`(c("deidentified_patient_id", 
                 "specimen_collection_date",
                 paste0(colnames(.)[3:ncol(.)], "_updated_Jan2025_YT"))) %>% 
  mutate(specimen_collection_date = as.Date(specimen_collection_date))

############################################################ III ### Replace old date and length/interval variable----
sequenced_patient_data <- sequenced_patient_data %>% 
  mutate(specimen_collection_date = as.Date(specimen_collection_date)) %>% 
  full_join(., updated_end_dates,
            by = c("deidentified_patient_id", 
                   "sequence_sample_vs_treatment" = "sequence_sample_vs_treatment_updated_Jan2025_YT", # I have to do that somehow otherwise there is a bug
                   "sample_treatment_sequence" = "sample_treatment_sequence_updated_Jan2025_YT",
                   "specimen_collection_date")) %>% 
  # add the variables at old location and check
  mutate(treatment_type_updated_newaddback = treatment_type_updated_Jan2025_YT, .after = treatment_type) %>% 
  mutate(date_of_first_corresponding_treatment_updated_newaddback = date_of_first_corresponding_treatment_updated_Jan2025_YT, .after = date_of_first_corresponding_treatment) %>% 
  mutate(chemotherapy_start_date_1_updated_newaddback = chemotherapy_start_date_1_updated_Jan2025_YT, .after = chemotherapy_start_date_1) %>% 
  mutate(chart_reviewed_chemotherapy_end_date_1_updated_updated_newaddback = chart_reviewed_chemotherapy_end_date_1_updated_updated_Jan2025_YT, .after = chart_reviewed_chemotherapy_end_date_1) %>% 
  mutate(interval_presample_chemo_updated_newaddback = interval_presample_chemo_updated_Jan2025_YT, .after = interval_presample_chemo) %>% 
  mutate(time_chemo_seqsample_updated_newaddback = time_chemo_seqsample_updated_Jan2025_YT, .after = time_chemo_seqsample) %>% 
  mutate(interval_prepost_sample_chemo_updated_newaddback = interval_prepost_sample_chemo_updated_Jan2025_YT, .after = interval_prepost_sample_chemo) %>% 
  mutate(chemo_length_indays_updated_newaddback = chemo_length_indays_updated_Jan2025_YT, .after = chemo_length_indays) %>% 
  mutate(radiation_start_date_1_updated_newaddback = radiation_start_date_1_updated_Jan2025_YT, .after = radiation_start_date_1) %>% 
  mutate(radiation_end_date1_updated_updated_newaddback = radiation_end_date1_updated_updated_Jan2025_YT, .after = radiation_end_date1_1) %>% 
  mutate(interval_presample_rad_updated_newaddback = interval_presample_rad_updated_Jan2025_YT, .after = interval_presample_rad) %>% 
  mutate(time_rad_seqsample_updated_newaddback = time_rad_seqsample_updated_Jan2025_YT, .after = time_rad_seqsample) %>% 
  mutate(interval_prepost_sample_rad_updated_newaddback = interval_prepost_sample_rad_updated_Jan2025_YT, .after = interval_prepost_sample_rad) %>% 
  mutate(radiation_length_indays_updated_newaddback = radiation_length_indays_updated_Jan2025_YT, .after = radiation_length_indays) %>% 
  mutate(hormone_therapy_start_date_1_updated_newaddback = hormone_therapy_start_date_1_updated_Jan2025_YT, .after = hormone_therapy_start_date_1) %>% 
  mutate(hormone_therapy_end_date1_1_updated_updated_newaddback = hormone_therapy_end_date1_1_updated_updated_Jan2025_YT, .after = hormone_therapy_end_date1_1) %>% 
  mutate(interval_presample_hormone_updated_newaddback = interval_presample_hormone_updated_Jan2025_YT, .after = interval_presample_hormone) %>% 
  mutate(time_hormone_seqsample_updated_newaddback = time_hormone_seqsample_updated_Jan2025_YT, .after = time_hormone_seqsample) %>% 
  mutate(interval_prepost_sample_hormone_updated_newaddback = interval_prepost_sample_hormone_updated_Jan2025_YT, .after = interval_prepost_sample_hormone) %>% 
  mutate(hormone_length_indays_updated_updated_newaddback = hormone_length_indays_updated_updated_Jan2025_YT, .after = hormone_length_indays) %>% 
  mutate(interval_presample_chemo_rad_updated_newaddback = interval_presample_chemo_rad_updated_Jan2025_YT, .after = interval_presample_chemo_rad) %>% 
  mutate(time_chemo_1st_seqsample_updated_newaddback = time_chemo_1st_seqsample_updated_Jan2025_YT, .after = time_chemo_1st_seqsample) %>% 
  mutate(time_1st_seqsample_rad_updated_newaddback = time_1st_seqsample_rad_updated_Jan2025_YT, .after = time_1st_seqsample_rad) %>% 
  mutate(time_rad_2nd_seqsample_updated_newaddback = time_rad_2nd_seqsample_updated_Jan2025_YT, .after = time_rad_2nd_seqsample) %>% 
  mutate(interval_prepost_sample_chemo_rad_updated_newaddback = interval_prepost_sample_chemo_rad_updated_Jan2025_YT, .after = interval_prepost_sample_chemo_rad) %>% 
  mutate(days_from_chemo_start_to_radiation_updated_newaddback = days_from_chemo_start_to_radiation_updated_Jan2025_YT, .after = days_from_chemo_start_to_radiation) %>% 
  mutate(chemotherapy_drug_1_updated_newaddback = chemotherapy_drug_1_updated_Jan2025_YT, .after = chemotherapy_drug_1) %>% 
  mutate(chemotherapy_drug_2_updated_newaddback = chemotherapy_drug_2_updated_Jan2025_YT, .after = chemotherapy_drug_2) %>% 
  mutate(chemotherapy_drug_3_updated_newaddback = chemotherapy_drug_3_updated_Jan2025_YT, .after = chemotherapy_drug_3) %>% 
  mutate(hormone_therapy_drug_1_updated_newaddback = hormone_therapy_drug_1_updated_Jan2025_YT, .after = hormone_therapy_drug_1) %>% 
  mutate(hormone_therapy_drug_2_updated_newaddback = hormone_therapy_drug_2_updated_Jan2025_YT, .after = hormone_therapy_drug_2) %>% 
  mutate(hormone_therapy_drug_3_updated_newaddback = hormone_therapy_drug_3_updated_Jan2025_YT, .after = hormone_therapy_drug_3) %>% 
  # Remove old and duplicate variables
  select(-ends_with("_updated_Jan2025_YT"),
         -c(treatment_type,
            date_of_first_corresponding_treatment, chemotherapy_start_date_1,
            radiation_start_date_1,
            chart_reviewed_chemotherapy_end_date_1, interval_presample_chemo,
            time_chemo_seqsample, interval_prepost_sample_chemo, chemo_length_indays,
            radiation_end_date1_1, interval_presample_rad,
            time_rad_seqsample, interval_prepost_sample_rad, radiation_length_indays,
            hormone_therapy_start_date_1, hormone_therapy_end_date1_1, interval_presample_hormone, 
            time_hormone_seqsample, interval_prepost_sample_hormone, 
            hormone_length_indays, interval_presample_chemo_rad, 
            time_chemo_1st_seqsample, time_1st_seqsample_rad, time_rad_2nd_seqsample,
            interval_prepost_sample_chemo_rad, days_from_chemo_start_to_radiation,
            chemotherapy_drug_1, chemotherapy_drug_2, chemotherapy_drug_3,
            hormone_therapy_drug_1, hormone_therapy_drug_2, hormone_therapy_drug_3)) %>%
  `colnames<-`(str_remove(colnames(.), "_updated_updated_newaddback|_updated_newaddback")) %>%
  # Update hormone date for a patient
  mutate(hormone_therapy_end_date1_1 = case_when(
    mrn == mrn_modif$mrn_hormone_end_date_is_na        ~ NA_Date_,
    TRUE                                               ~ hormone_therapy_end_date1_1
  )) %>% 
  # Update the combined end date
  mutate(chart_reviewed_end_date_of_first_corresponding_treatment = case_when(
    treatment_type == "chemo"              ~ chart_reviewed_chemotherapy_end_date_1,
    treatment_type == "chemorad"           ~ chart_reviewed_chemotherapy_end_date_1,
    treatment_type == "hormone"            ~ hormone_therapy_end_date1_1,
    treatment_type == "radiation"          ~ radiation_end_date1
  ))

  
############################################################ IV ### Code NADIR and Cytopenia----
# Create a data with date range
temp <- sequenced_patient_data %>% 
  select(mrn, radiation_end_date1, 
         hormone_therapy_end_date1_1,
         chart_reviewed_chemotherapy_end_date_1,
         treatment_type
  ) %>% 
  distinct() %>% 
  mutate(treatment_enddate = case_when(
    treatment_type == "radiation"                            ~ radiation_end_date1,
    treatment_type == "hormone"                              ~ hormone_therapy_end_date1_1,
    str_detect(treatment_type, "chemo")                      ~ chart_reviewed_chemotherapy_end_date_1
  )) %>% 
  mutate(max_persistent_cytopenia_date = case_when(
    treatment_type == "radiation"                            
    ~ add_with_rollback(treatment_enddate , months(1), roll_to_first = TRUE),
    
    treatment_type == "hormone"                  
    ~ add_with_rollback(treatment_enddate , months(1), roll_to_first = TRUE),
    
    str_detect(treatment_type, "chemo")                      
    ~ add_with_rollback(treatment_enddate , months(1), roll_to_first = TRUE)
    
  )) %>% 
  mutate(max_prolonged_cytopenia_date = case_when(
    treatment_type == "radiation"                            
    ~ add_with_rollback(treatment_enddate , months(3), roll_to_first = TRUE),
    
    treatment_type == "hormone"                 
    ~ add_with_rollback(treatment_enddate , months(3), roll_to_first = TRUE),
    
    str_detect(treatment_type, "chemo")                     
    ~ add_with_rollback(treatment_enddate , months(3), roll_to_first = TRUE)
    
  )) %>% 
  mutate(reason_for_no_persistent_prolonged_cytopenia_date = case_when(
    treatment_type == "radiation" &
      is.na(radiation_end_date1)                                    ~ "no radiation end date",
    treatment_type == "hormone" &
      is.na(hormone_therapy_end_date1_1)                            ~ "no hormone end date",
    str_detect(treatment_type, "chemo") &
      is.na(chart_reviewed_chemotherapy_end_date_1)                 ~ "no chemo end date"
  ))

# Cytopenia----
# Cytopenias were defined by using World Health Organization criteria
# anemia = hemoglobin concentration <13.0 g/dl in male participants 
# and <12.0 g/dl in female participants; 
# thrombocytopenia = platelet counts <150 × 109 cells/l; 
# and neutropenia = absolute neutrophil count <1.8 × 109 cells/l

# cells/ l == cells / 1000ul == cells /ul
# 
# 1000cells /ul == 1000*1000cells / 1000ul == 1000+1000cells / L == 1 x 106


# HGB
hgb <- hgb_plt %>%
  mutate(mrn = as.character(mrn)) %>% 
  select(mrn, lab_nm, lab_result, lab_unit, lab_date) %>% 
  filter(lab_nm == "Hemoglobin") %>% 
  mutate(lab_date = as.POSIXct(lab_date)) %>% 
  right_join(., temp, by = "mrn") %>% 
  mutate(persistent_low_hgb = case_when(
    is.na(max_persistent_cytopenia_date)                    ~ "No treatment end date",
    lab_date >= treatment_enddate &
      lab_date <= max_persistent_cytopenia_date &
      lab_result >= 12                                      ~ "No" # anemia is <12 is for female
  )) %>% 
  mutate(prolonged_low_hgb = case_when(
    is.na(max_prolonged_cytopenia_date)                     ~ "No treatment end date",
    lab_date >= treatment_enddate &
      lab_date <= max_prolonged_cytopenia_date &
      lab_result >= 12                                      ~ "No"
  )) %>% 
  mutate(persistent_low_hgb_buffertime = case_when(
    is.na(max_persistent_cytopenia_date)                    ~ "No treatment end date",
    lab_date >= treatment_enddate &
      lab_date <= max_persistent_cytopenia_date + days(14) &
      lab_date > max_persistent_cytopenia_date &
      lab_result >= 12                                      ~ "No" # anemia is <12 is for female
  )) %>% 
  mutate(prolonged_low_hgb_buffertime = case_when(
    is.na(max_prolonged_cytopenia_date)                     ~ "No treatment end date",
    lab_date >= treatment_enddate &
      lab_date <= max_prolonged_cytopenia_date + days(14) &
      lab_date > max_prolonged_cytopenia_date &
      lab_result >= 12                                      ~ "No"
  )) %>% 
  group_by(mrn) %>%
  fill(persistent_low_hgb, prolonged_low_hgb,
       persistent_low_hgb_buffertime,
       prolonged_low_hgb_buffertime, .direction = "updown") %>%
  ungroup() %>%
  mutate(persistent_low_hgb = case_when(
    persistent_low_hgb == "No"                              ~ "No",
    persistent_low_hgb == "No treatment end date"           ~ "No treatment end date",
    lab_date < treatment_enddate |
      lab_date > max_persistent_cytopenia_date              ~ "Not correct lab date",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(persistent_low_hgb)                               ~ "Yes"
  )) %>% 
  mutate(prolonged_low_hgb = case_when(
    prolonged_low_hgb == "No"                               ~ "No",
    prolonged_low_hgb == "No treatment end date"            ~ "No treatment end date",
    lab_date < treatment_enddate |
      lab_date > max_prolonged_cytopenia_date               ~ "Not correct lab date",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(prolonged_low_hgb)                                ~ "Yes"
  )) %>% 
  mutate(persistent_low_hgb_buffertime = case_when(
    persistent_low_hgb_buffertime == "No"                              ~ "No",
    persistent_low_hgb_buffertime == "No treatment end date"           ~ "No treatment end date",
    lab_date <= max_persistent_cytopenia_date + days(14) |
      lab_date > max_persistent_cytopenia_date              ~ "Not correct lab date",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(persistent_low_hgb_buffertime)                               ~ "Yes"
  )) %>% 
  mutate(prolonged_low_hgb_buffertime = case_when(
    prolonged_low_hgb_buffertime == "No"                               ~ "No",
    prolonged_low_hgb_buffertime == "No treatment end date"            ~ "No treatment end date",
    lab_date <= max_persistent_cytopenia_date + days(14) |
      lab_date > max_persistent_cytopenia_date               ~ "Not correct lab date",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(prolonged_low_hgb_buffertime)                                ~ "Yes"
  )) %>% 
  select(mrn, persistent_low_hgb, prolonged_low_hgb, 
         persistent_low_hgb_buffertime, prolonged_low_hgb_buffertime) %>% 
  mutate(never_correct_lab_date = case_when(
    persistent_low_hgb_buffertime == "Not correct lab date" & 
      prolonged_low_hgb_buffertime == "Not correct lab date" &
      persistent_low_hgb == "Not correct lab date" & 
      prolonged_low_hgb == "Not correct lab date"               ~ "Never correct lab date"
  )) %>% 
  mutate_at(c("persistent_low_hgb",
              "prolonged_low_hgb",
              "persistent_low_hgb_buffertime",
              "prolonged_low_hgb_buffertime"), ~ na_if(., "Not correct lab date")) %>%
  group_by(mrn) %>%
    fill(persistent_low_hgb, prolonged_low_hgb,
         persistent_low_hgb_buffertime,
         prolonged_low_hgb_buffertime, .direction = "updown") %>%
  ungroup() %>%
  distinct(mrn, .keep_all = TRUE) %>% 
  mutate(persistent_low_hgb = coalesce(persistent_low_hgb, persistent_low_hgb_buffertime, never_correct_lab_date)) %>% 
  mutate(prolonged_low_hgb = coalesce(prolonged_low_hgb, prolonged_low_hgb_buffertime, never_correct_lab_date)) %>% 
  select(mrn, persistent_low_hgb, prolonged_low_hgb)

# PLT
plt <- hgb_plt %>%
  mutate(mrn = as.character(mrn)) %>% 
  select(mrn, lab_nm, lab_result, lab_unit, lab_date) %>% 
  filter(lab_nm == "Platelet Count(k/uL)") %>% 
  mutate(lab_date = as.POSIXct(lab_date)) %>% 
  right_join(., temp, by = "mrn") %>% 
  mutate(persistent_low_plt = case_when(
    is.na(max_persistent_cytopenia_date)                    ~ "No treatment end date",
    lab_date >= treatment_enddate &
      lab_date <= max_persistent_cytopenia_date &
      lab_result >= 150                                      ~ "No" # thrombocytopenia = plt <150 × 109 cells/l
  )) %>% 
  mutate(prolonged_low_plt = case_when(
    is.na(max_prolonged_cytopenia_date)                     ~ "No treatment end date",
    lab_date >= treatment_enddate &
      lab_date <= max_prolonged_cytopenia_date &
      lab_result >= 150                                      ~ "No"
  )) %>% 
  mutate(persistent_low_plt_buffertime = case_when(
    is.na(max_persistent_cytopenia_date)                    ~ "No treatment end date",
    lab_date >= treatment_enddate &
      lab_date <= max_persistent_cytopenia_date + days(14) &
      lab_date > max_persistent_cytopenia_date &
      lab_result >= 12                                      ~ "No" # anemia is <12 is for female
  )) %>% 
  mutate(prolonged_low_plt_buffertime = case_when(
    is.na(max_prolonged_cytopenia_date)                     ~ "No treatment end date",
    lab_date >= treatment_enddate &
      lab_date <= max_prolonged_cytopenia_date + days(14) &
      lab_date > max_prolonged_cytopenia_date &
      lab_result >= 12                                      ~ "No"
  )) %>% 
  group_by(mrn) %>%
  fill(persistent_low_plt, prolonged_low_plt,
       persistent_low_plt_buffertime,
       prolonged_low_plt_buffertime, .direction = "updown") %>%
  ungroup() %>%
  mutate(persistent_low_plt = case_when(
    persistent_low_plt == "No"                              ~ "No",
    persistent_low_plt == "No treatment end date"           ~ "No treatment end date",
    lab_date < treatment_enddate |
      lab_date > max_persistent_cytopenia_date              ~ "Not correct lab date",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(persistent_low_plt)                               ~ "Yes"
  )) %>% 
  mutate(prolonged_low_plt = case_when(
    prolonged_low_plt == "No"                               ~ "No",
    prolonged_low_plt == "No treatment end date"            ~ "No treatment end date",
    lab_date < treatment_enddate |
      lab_date > max_prolonged_cytopenia_date               ~ "Not correct lab date",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(prolonged_low_plt)                                ~ "Yes"
  )) %>% 
  mutate(persistent_low_plt_buffertime = case_when(
    persistent_low_plt_buffertime == "No"                              ~ "No",
    persistent_low_plt_buffertime == "No treatment end date"           ~ "No treatment end date",
    lab_date <= max_persistent_cytopenia_date + days(14) |
      lab_date > max_persistent_cytopenia_date              ~ "Not correct lab date",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(persistent_low_plt_buffertime)                               ~ "Yes"
  )) %>% 
  mutate(prolonged_low_plt_buffertime = case_when(
    prolonged_low_plt_buffertime == "No"                               ~ "No",
    prolonged_low_plt_buffertime == "No treatment end date"            ~ "No treatment end date",
    lab_date <= max_persistent_cytopenia_date + days(14) |
      lab_date > max_persistent_cytopenia_date               ~ "Not correct lab date",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(prolonged_low_plt_buffertime)                                ~ "Yes"
  )) %>% 
  select(mrn, persistent_low_plt, prolonged_low_plt, 
         persistent_low_plt_buffertime, prolonged_low_plt_buffertime) %>% 
  mutate(never_correct_lab_date = case_when(
    persistent_low_plt_buffertime == "Not correct lab date" & 
      prolonged_low_plt_buffertime == "Not correct lab date" &
      persistent_low_plt == "Not correct lab date" & 
      prolonged_low_plt == "Not correct lab date"               ~ "Never correct lab date"
  )) %>% 
  mutate_at(c("persistent_low_plt",
              "prolonged_low_plt",
              "persistent_low_plt_buffertime",
              "prolonged_low_plt_buffertime"), ~ na_if(., "Not correct lab date")) %>%
  group_by(mrn) %>%
  fill(persistent_low_plt, prolonged_low_plt,
       persistent_low_plt_buffertime,
       prolonged_low_plt_buffertime, .direction = "updown") %>%
  ungroup() %>%
  distinct(mrn, .keep_all = TRUE) %>% 
  mutate(persistent_low_plt = coalesce(persistent_low_plt, persistent_low_plt_buffertime, never_correct_lab_date)) %>% 
  mutate(prolonged_low_plt = coalesce(prolonged_low_plt, prolonged_low_plt_buffertime, never_correct_lab_date)) %>% 
  select(mrn, persistent_low_plt, prolonged_low_plt)

# ANC
# I take neutro or auto first and when missing add the others...
neutrophil_other <- neutrophil %>% 
  mutate(mrn = as.character(mrn)) %>% 
  # select neutroplil poly and bands
  filter(lab_nm != "Neutrophil") %>% 
  select(mrn, lab_nm, lab_result, lab_unit, lab_neutrophil_date = order_dtm) %>% 
  mutate(lab_result = case_when(
    lab_unit == "k/uL"           ~ lab_result,
    lab_unit == "cells/uL"       ~ lab_result / 1000
  )) %>% 
  # Multiple values on the same day - pick the lowest
  arrange(mrn, lab_neutrophil_date, lab_nm) %>% 
  distinct(mrn, lab_neutrophil_date, lab_nm, .keep_all = TRUE) %>% 
  # pivot and add ploy + bands values together
  mutate(lab_nm = str_replace(lab_nm, " ", "_"),
         lab_nm = str_to_lower(lab_nm)) %>% 
  pivot_wider(id_cols = c(mrn, lab_neutrophil_date), 
              names_from = lab_nm, 
              values_from = c(lab_result, lab_unit), 
              names_vary = "slowest") %>% 
  mutate(lab_result_neutrophil = lab_result_neutrophil_bands + lab_result_neutrophil_poly) %>% 
  filter(!is.na(lab_result_neutrophil))

neutrophil <- neutrophil %>% 
  mutate(mrn = as.character(mrn)) %>% 
  filter(lab_nm == "Neutrophil") %>% 
  select(mrn, lab_result_neutrophil = lab_result, 
         lab_neutrophil_unit = lab_unit, 
         lab_neutrophil_date = order_dtm) %>% 
  bind_rows(., neutrophil_other) %>% 
  # I verify all units of selected values are k/uL
  mutate(lab_neutrophil_unit = "k/uL")

anc <- neutrophil %>%
  mutate(lab_neutrophil_date = as.POSIXct(lab_neutrophil_date)) %>% 
  right_join(., temp, by = "mrn") %>% 
  mutate(persistent_low_anc = case_when(
    is.na(max_persistent_cytopenia_date)                    ~ "No treatment end date",
    lab_neutrophil_date >= treatment_enddate &
      lab_neutrophil_date <= max_persistent_cytopenia_date &
      lab_result_neutrophil >= 1.8                          ~ "No" # neutropenia = anc <1.8 × 109 cells/l
  )) %>% 
  mutate(prolonged_low_anc = case_when(
    is.na(max_prolonged_cytopenia_date)                     ~ "No treatment end date",
    lab_neutrophil_date >= treatment_enddate &
      lab_neutrophil_date <= max_prolonged_cytopenia_date &
      lab_result_neutrophil >= 1.8                          ~ "No"
  )) %>% 
  mutate(persistent_low_anc_buffertime = case_when(
    is.na(max_persistent_cytopenia_date)                    ~ "No treatment end date",
    lab_neutrophil_date >= treatment_enddate &
      lab_neutrophil_date <= max_persistent_cytopenia_date + days(14) &
      lab_neutrophil_date > max_persistent_cytopenia_date &
      lab_result_neutrophil >= 12                                      ~ "No" # anemia is <12 is for female
  )) %>% 
  mutate(prolonged_low_anc_buffertime = case_when(
    is.na(max_prolonged_cytopenia_date)                     ~ "No treatment end date",
    lab_neutrophil_date >= treatment_enddate &
      lab_neutrophil_date <= max_prolonged_cytopenia_date + days(14) &
      lab_neutrophil_date > max_prolonged_cytopenia_date &
      lab_result_neutrophil >= 12                                      ~ "No"
  )) %>% 
  group_by(mrn) %>%
  fill(persistent_low_anc, prolonged_low_anc,
       persistent_low_anc_buffertime,
       prolonged_low_anc_buffertime, .direction = "updown") %>%
  ungroup() %>%
  mutate(persistent_low_anc = case_when(
    persistent_low_anc == "No"                              ~ "No",
    persistent_low_anc == "No treatment end date"           ~ "No treatment end date",
    lab_neutrophil_date < treatment_enddate |
      lab_neutrophil_date > max_persistent_cytopenia_date   ~ "Not correct lab date",
    is.na(lab_result_neutrophil)                            ~ NA_character_,
    is.na(persistent_low_anc)                               ~ "Yes"
  )) %>% 
  mutate(prolonged_low_anc = case_when(
    prolonged_low_anc == "No"                               ~ "No",
    prolonged_low_anc == "No treatment end date"            ~ "No treatment end date",
    lab_neutrophil_date < treatment_enddate |
      lab_neutrophil_date > max_prolonged_cytopenia_date    ~ "Not correct lab date",
    is.na(lab_result_neutrophil)                            ~ NA_character_,
    is.na(prolonged_low_anc)                                ~ "Yes"
  )) %>% 
  mutate(persistent_low_anc_buffertime = case_when(
    persistent_low_anc_buffertime == "No"                              ~ "No",
    persistent_low_anc_buffertime == "No treatment end date"           ~ "No treatment end date",
    lab_neutrophil_date <= max_persistent_cytopenia_date + days(14) |
      lab_neutrophil_date > max_persistent_cytopenia_date              ~ "Not correct lab date",
    is.na(lab_result_neutrophil)                                       ~ NA_character_,
    is.na(persistent_low_anc_buffertime)                               ~ "Yes"
  )) %>% 
  mutate(prolonged_low_anc_buffertime = case_when(
    prolonged_low_anc_buffertime == "No"                               ~ "No",
    prolonged_low_anc_buffertime == "No treatment end date"            ~ "No treatment end date",
    lab_neutrophil_date <= max_persistent_cytopenia_date + days(14) |
      lab_neutrophil_date > max_persistent_cytopenia_date               ~ "Not correct lab date",
    is.na(lab_result_neutrophil)                                       ~ NA_character_,
    is.na(prolonged_low_anc_buffertime)                                ~ "Yes"
  )) %>% 
  select(mrn, persistent_low_anc, prolonged_low_anc, 
         persistent_low_anc_buffertime, prolonged_low_anc_buffertime) %>% 
  mutate(never_correct_lab_date = case_when(
    persistent_low_anc_buffertime == "Not correct lab date" & 
      prolonged_low_anc_buffertime == "Not correct lab date" &
      persistent_low_anc == "Not correct lab date" & 
      prolonged_low_anc == "Not correct lab date"               ~ "Never correct lab date"
  )) %>% 
  mutate_at(c("persistent_low_anc",
              "prolonged_low_anc",
              "persistent_low_anc_buffertime",
              "prolonged_low_anc_buffertime"), ~ na_if(., "Not correct lab date")) %>%
  group_by(mrn) %>%
  fill(persistent_low_anc, prolonged_low_anc,
       persistent_low_anc_buffertime,
       prolonged_low_anc_buffertime, .direction = "updown") %>%
  ungroup() %>%
  distinct(mrn, .keep_all = TRUE) %>% 
  mutate(persistent_low_anc = coalesce(persistent_low_anc, persistent_low_anc_buffertime, never_correct_lab_date)) %>% 
  mutate(prolonged_low_anc = coalesce(prolonged_low_anc, prolonged_low_anc_buffertime, never_correct_lab_date)) %>% 
  select(mrn, persistent_low_anc, prolonged_low_anc)

sequenced_patient_data <- sequenced_patient_data %>% 
  left_join(., hgb, by = "mrn") %>% 
  left_join(., plt, by = "mrn") %>% 
  left_join(., anc, by = "mrn") %>% 
  mutate(persistent_cytopenia = case_when(
    persistent_low_hgb == "Yes" |
      persistent_low_plt == "Yes" |
      persistent_low_anc == "Yes"                            ~ "Yes",
    persistent_low_hgb == "No" |
      persistent_low_plt == "No" |
      persistent_low_anc == "No"                             ~ "No",
  )) %>% 
  mutate(prolonged_cytopenia = case_when(
    prolonged_low_hgb == "Yes" |
      prolonged_low_plt == "Yes" |
      prolonged_low_anc == "Yes"                             ~ "Yes",
    prolonged_low_hgb == "No" |
      prolonged_low_plt == "No" |
      prolonged_low_anc == "No"                              ~ "No",
  ))

write_csv(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia_",
                                         today(), ".csv"))
write_csv(sequenced_patient_data, paste0(path_save, "/processed data",
                                         "/Identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia_",
                                         today(), ".csv"))


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
                 "/De-identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia_", 
                 today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data", 
                 "/De-identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia_", 
                 today(), ".csv"))


# clean up
rm(temp, hgb, plt, anc,
   neutrophil_other, updated_end_dates,
   deids_data)

# NADIR----
temp <- sequenced_patient_data %>% 
  select(mrn, 
         date_of_first_corresponding_treatment) %>% 
  distinct() %>% 
  mutate(date_of_first_corresponding_treatment = as.POSIXct(date_of_first_corresponding_treatment))

chart_dat <- chart_dat %>% 
  mutate(mrn = as.character(mrn)) %>% 
  select(mrn, 
         starts_with("wbc_nadir"),
         starts_with("rbc_nadir"),
         starts_with("hgb_nadir")) %>% 
  inner_join(., temp, 
             by = "mrn") %>% 
  mutate(keep = case_when(
    wbc_nadir_date >= (date_of_first_corresponding_treatment + days(7)) &
      wbc_nadir_date <= (date_of_first_corresponding_treatment + days(15))           ~ "Yes"
  )) %>% 
  mutate_at(c("wbc_nadir",
              "wbc_nadir_date"), ~ case_when(
                keep == "Yes"           ~ .,
                is.na(keep)             ~ NA)) %>% 
  mutate(keep = case_when(
    rbc_nadir_date >= (date_of_first_corresponding_treatment + days(7)) &
      rbc_nadir_date <= (date_of_first_corresponding_treatment + days(15))           ~ "Yes"
  )) %>% 
  mutate_at(c("rbc_nadir",
              "rbc_nadir_date"), ~ case_when(
                keep == "Yes"           ~ .,
                is.na(keep)             ~ NA)) %>% 
  mutate(keep = case_when(
    hgb_nadir_date >= (date_of_first_corresponding_treatment + days(7)) &
      hgb_nadir_date <= (date_of_first_corresponding_treatment + days(15))           ~ "Yes"
  )) %>% 
  mutate_at(c("hgb_nadir",
              "hgb_nadir_date"), ~ case_when(
                keep == "Yes"           ~ .,
                is.na(keep)             ~ NA)) %>% 
  filter(!is.na(wbc_nadir) | 
           !is.na(rbc_nadir) |
           !is.na(hgb_nadir)) %>% 
  select(mrn, 
         starts_with("wbc_nadir"),
         starts_with("rbc_nadir"),
         starts_with("hgb_nadir"))

cbc1 <- cbc %>% 
  select(patient_id, lab_nm, lab_result, lab_unit, lab_date) %>% 
  mutate(lab_result = as.numeric(lab_result)) %>% 
  #add mrn
  left_join(., Demographic %>% 
              select(mrn, patient_id), by = c("patient_id")) %>% 
  filter((lab_nm == "Hemoglobin" & lab_unit == "g/dL") | 
           (lab_nm == "Platelet Count(k/uL)" & lab_unit == "k/uL") |
           (lab_nm == "WBC(k/uL)" & lab_unit == "k/uL") |
           lab_nm == "MCV" | lab_nm == "RDW" | lab_nm == "RBC") %>% 
  mutate(lab_nm = case_when(
    lab_nm == "Hemoglobin"            ~ "hgb",
    lab_nm == "Platelet Count(k/uL)"  ~ "plt",
    lab_nm == "WBC(k/uL)"             ~ "wbc",
    lab_nm == "MCV"                   ~ "mcv",
    lab_nm == "RDW"                   ~ "rdw",
    lab_nm == "RBC"                   ~ "rbc"
  )) %>% 
  # add date of pre-sample and select closest cbc
  inner_join(., temp, 
             by = "mrn") %>% 
  # filter date before treatment
  mutate(keep = case_when(
    lab_date >= (date_of_first_corresponding_treatment + days(7)) &
      lab_date <= (date_of_first_corresponding_treatment + days(15))           ~ "Yes"
  )) %>% 
  filter(keep == "Yes") %>% 
  # arrange to get the lowest if multiple
  arrange(mrn, lab_nm, lab_result) %>% 
  distinct(mrn, lab_nm, .keep_all = TRUE) %>% 
  # pivot
  select(mrn, lab_nm, lab_result, lab_unit, lab_date) %>% 
  pivot_wider(id_cols = mrn, 
              names_from = lab_nm, 
              values_from = c(lab_result, lab_unit, lab_date), 
              names_vary = "slowest") %>% 
  `colnames<-`(c("mrn", 
                 paste0("nadir_", colnames(.)[2:ncol(.)]))) %>% 
  full_join(., chart_dat, by = "mrn") %>% 
  inner_join(., temp, 
             by = "mrn") %>% 
  mutate(nadir_lab_result_hgb = coalesce(hgb_nadir, nadir_lab_result_hgb)) %>% 
  mutate(nadir_lab_date_hgb = coalesce(hgb_nadir_date, nadir_lab_date_hgb)) %>% 
  mutate(nadir_lab_result_rbc = coalesce(rbc_nadir, nadir_lab_result_rbc)) %>% 
  mutate(nadir_lab_date_rbc = coalesce(rbc_nadir_date, nadir_lab_date_rbc)) %>% 
  mutate(nadir_lab_result_wbc = coalesce(wbc_nadir, nadir_lab_result_wbc)) %>% 
  mutate(nadir_lab_date_wbc = coalesce(wbc_nadir_date, nadir_lab_date_wbc)) %>% 
  select(mrn, starts_with("nadir_"))

neutrophil <- neutrophil %>%
  select(mrn, lab_result_neutrophil, lab_neutrophil_unit, lab_neutrophil_date) %>% 
  right_join(., temp, by = "mrn") %>% 
  mutate(keep = case_when(
    lab_neutrophil_date >= (date_of_first_corresponding_treatment + days(7)) &
      lab_neutrophil_date <= (date_of_first_corresponding_treatment + days(15))           ~ "Yes"
  )) %>% 
  filter(keep == "Yes") %>% 
  # arrange to get the lowest if multiple
  arrange(mrn, lab_result_neutrophil) %>% 
  # pick the lowest
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, nadir_lab_result_anc = lab_result_neutrophil, 
         nadir_lab_unit_anc = lab_neutrophil_unit, 
         nadir_lab_date_anc = lab_neutrophil_date)

sequenced_patient_data <- sequenced_patient_data %>% 
  left_join(., cbc1, by = "mrn") %>% 
  left_join(., neutrophil, by = "mrn")

write_csv(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia-NADIR_",
                                         today(), ".csv"))
write_csv(sequenced_patient_data, paste0(path_save, "/processed data",
                                         "/Identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia-NADIR_",
                                         today(), ".csv"))


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
                 "/De-identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia-NADIR_", 
                 today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data", 
                 "/De-identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia-NADIR_", 
                 today(), ".csv"))










