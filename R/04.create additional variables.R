# Import library
library(tidyverse)

# Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                 "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read.csv(paste0(#path,
  here::here(),
  "/processed data/Identified breast data with sequenced sequential sample and clinical_07092024.csv"))
wbc <- read.csv(paste0(#path,
  here::here(),
  "/processed data/Cleaned WBC data.csv"))
hgb_plt <- read.csv(paste0(#path,
  here::here(),
  "/processed data/Cleaned HgB-Plt data.csv"))
neutrophil <- 
  read_csv(paste0(path, "/processed data/Neutrophil lab data.csv"))

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


# CBC
head(cbc)
table(cbc$lab_nm)
Have RBC, hemoglobin, and hematocrit                                  HCT                           Hemoglobin                              MCH(pg) 
111                                50266                                50269                                50156 
MCHC                                  MCV                                  MPV                 Platelet Count(k/uL) 
50156                                50156                                31106                                50314 
RBC                                  RDW                            WBC(k/uL) 
50156                                49880                                50156 
dat <- sequenced_patient_data %>% 
  mutate(mrn = as.character(mrn)) %>% 
  mutate(presample_date = case_when(
    sample_treatment_sequence == "pre"                       ~ specimen_collection_date
  )) %>% 
  mutate(postsample_date = case_when(
    sample_treatment_sequence == "post"                      ~ specimen_collection_date
  )) %>% 
  group_by(mrn) %>% 
  fill(presample_date, postsample_date, .direction = "updown") %>% 
  ungroup() %>% 
  select(mrn, presample_date, postsample_date) %>% 
  distinct(mrn, .keep_all = TRUE)

RBC <- cbc %>% 
  filter(lab_nm == "RBC") %>% 
  mutate(lab_result = as.numeric(lab_result)) %>% 
  mutate(lab_result = case_when(
    lab_unit == "mil/uL"                                     ~ lab_result,
    # lab_unit == "mil/mL"                                     ~ lab_result / 1000
  )) %>% 
  mutate(lab_unit = "mil/uL") %>% 
  #add mrn
  left_join(., Demographic %>% 
              select(mrn, patient_id), by = c("patient_id")) %>% 
  select(mrn, lab_result_rbc = lab_result, lab_unit_rbc = lab_unit, lab_date_rbc = order_dtm) # NEED TO ASK UNITS -------------
rbc <- RBC %>% 
  # Add date's range
  right_join(., dat, by = "mrn") %>% 
  mutate(keep = case_when(
    lab_date_rbc >= presample_date &
      lab_date_rbc <= postsample_date                        ~ "Yes"
  )) %>% 
  filter(keep == "Yes") %>% 
  arrange(mrn, lab_result_rbc) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, lowest_rbc_lab_value = lab_result_rbc, 
         lowest_rbc_lab_unit = lab_unit_rbc, 
         lowest_rbc_lab_date = lab_date_rbc)

mcv <- cbc %>% 
  filter(lab_nm == "MCV") %>% 
  mutate(lab_result = as.numeric(lab_result)) %>% 
  #add mrn
  left_join(., Demographic %>% 
              select(mrn, patient_id), by = c("patient_id")) %>% 
  select(mrn, lab_result_mcv = lab_result, lab_unit_mcv = lab_unit, lab_date_mcv = order_dtm) %>% 
  # Add date's range
  right_join(., dat, by = "mrn") %>% 
  mutate(keep = case_when(
    lab_date_mcv >= presample_date &
      lab_date_mcv <= postsample_date                        ~ "Yes"
  )) %>% 
  filter(keep == "Yes") %>% 
  arrange(mrn, lab_result_mcv) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, lowest_mcv_lab_value = lab_result_mcv, 
         lowest_mcv_lab_unit = lab_unit_mcv, 
         lowest_mcv_lab_date = lab_date_mcv)

rdw <- cbc %>% 
  filter(lab_nm == "RDW") %>% 
  mutate(lab_result = as.numeric(lab_result)) %>% 
  #add mrn
  left_join(., Demographic %>% 
              select(mrn, patient_id), by = c("patient_id")) %>% 
  select(mrn, lab_result_rdw = lab_result, lab_unit_rdw = lab_unit, lab_date_rdw = order_dtm) %>% 
  # Add date's range
  right_join(., dat, by = "mrn") %>% 
  mutate(keep = case_when(
    lab_date_rdw >= presample_date &
      lab_date_rdw <= postsample_date                        ~ "Yes"
  )) %>% 
  filter(keep == "Yes") %>% 
  arrange(mrn, lab_result_rdw) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, lowest_rdw_lab_value = lab_result_rdw, 
         lowest_rdw_lab_unit = lab_unit_rdw, 
         lowest_rdw_lab_date = lab_date_rdw)

wbc <- wbc %>% 
  mutate(mrn = as.character(mrn)) %>% 
  # Add date's range
  right_join(., dat, by = "mrn") %>% 
  mutate(keep = case_when(
    lab_date_wbc >= presample_date &
      lab_date_wbc <= postsample_date                        ~ "Yes"
  )) %>% 
  filter(keep == "Yes") %>% 
  arrange(mrn, lab_result_wbc) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, lowest_wbc_lab_value = lab_result_wbc, 
         lowest_wbc_lab_unit = lab_unit_wbc, 
         lowest_wbc_lab_date = lab_date_wbc) %>% 
  mutate(lowest_wbc_lab_date = as.POSIXct(lowest_wbc_lab_date))

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
  select(mrn, lab_result_neutrophil, lab_neutrophil_unit, lab_neutrophil_date) %>% 
  # Add date's range
  right_join(., dat, by = "mrn") %>% 
  mutate(keep = case_when(
    lab_neutrophil_date >= presample_date &
      lab_neutrophil_date <= postsample_date                 ~ "Yes"
  )) %>% 
  filter(keep == "Yes") %>% 
  arrange(mrn, lab_result_neutrophil) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, lowest_anc_lab_value = lab_result_neutrophil, 
         lowest_anc_lab_unit = lab_neutrophil_unit, 
         lowest_anc_lab_date = lab_neutrophil_date)

HGB_PLT <- hgb_plt %>% 
  mutate(mrn = as.character(mrn)) %>% 
  select(mrn, lab_nm, lab_result, lab_unit, lab_date)
hgb_plt <- HGB_PLT %>% 
  # Add date's range
  right_join(., dat, by = "mrn") %>% 
  mutate(keep = case_when(
    lab_date >= presample_date &
      lab_date <= postsample_date                            ~ "Yes"
  )) %>% 
  filter(keep == "Yes") %>% 
  arrange(mrn, lab_nm, lab_result) %>% 
  distinct(mrn, lab_nm, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = mrn, 
              names_from = lab_nm, 
              values_from = c(lab_result, lab_unit, lab_date), 
              names_vary = "slowest") %>% 
  select(mrn, lowest_hgb_lab_value = lab_result_Hemoglobin, 
         lowest_hgb_lab_unit = lab_unit_Hemoglobin, 
         lowest_hgb_lab_date = lab_date_Hemoglobin,
         lowest_plt_lab_value = `lab_result_Platelet Count(k/uL)`, 
         lowest_plt_lab_unit = `lab_unit_Platelet Count(k/uL)`, 
         lowest_plt_lab_date = `lab_date_Platelet Count(k/uL)`) %>% 
  mutate(lowest_hgb_lab_date = as.POSIXct(lowest_hgb_lab_date),
         lowest_plt_lab_date = as.POSIXct(lowest_plt_lab_date))
  

lowest_cbc <- full_join(rbc, wbc, by = "mrn") %>% 
  full_join(., anc, by = "mrn") %>% 
  full_join(., hgb_plt, by = "mrn") %>% 
  full_join(., mcv, by = "mrn") %>% 
  full_join(., rdw, by = "mrn")
  

sequenced_patient_data <- sequenced_patient_data %>% 
  mutate(mrn = as.character(mrn)) %>% 
  left_join(., lowest_cbc, by = "mrn")

write_csv(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data with sequenced sequential sample and clinical-CBC_",
                                         today(), ".csv"))
write_csv(sequenced_patient_data, paste0(path, "/processed data",
                                         "/Identified breast data with sequenced sequential sample and clinical-CBC_",
                                         today(), ".csv"))

# Cytopenia
# Cytopenias were defined by using World Health Organization criteria ######### Do it later
# anemia = hemoglobin concentration <13.0 g/dl in male participants 
# and <12.0 g/dl in female participants; 
# thrombocytopenia = platelet counts <150 × 109 cells/l; 
# and neutropenia = absolute neutrophil count <1.8 × 109 cells/l

# cells/ l == cells / 1000ul == cells /ul
# 
# 1000cells /ul == 1000*1000cells / 1000ul == 1000+1000cells / L == 1 x 106

temp <- sequenced_patient_data %>% 
  select(mrn, date_of_first_corresponding_treatment, 
         hormone_therapy_end_date1_1,
         chemotherapy_end_date1_1,
         treatment_type
         ) %>% 
  distinct() %>% 
  mutate_at(c("date_of_first_corresponding_treatment",
              "hormone_therapy_end_date1_1",
              "chemotherapy_end_date1_1"), ~ as.POSIXct(.)) %>% 
  mutate(max_persistent_cytopenia_date = case_when(
    treatment_type == "radiation"                            ~ date_of_first_corresponding_treatment + months(1),
    
    treatment_type == "hormone" &
      is.na(hormone_therapy_end_date1_1)                     ~ NA_POSIXct_,
    treatment_type == "hormone" &
      !is.na(hormone_therapy_end_date1_1)                    ~ hormone_therapy_end_date1_1 + months(1),
    
    str_detect(treatment_type, "chemo") &
      is.na(chemotherapy_end_date1_1)                        ~ NA_POSIXct_,
    str_detect(treatment_type, "chemo") &
      !is.na(chemotherapy_end_date1_1)                       ~ chemotherapy_end_date1_1 + months(1),
    
  )) %>% 
  mutate(max_prolonged_cytopenia_date = case_when(
    treatment_type == "radiation"                            ~ date_of_first_corresponding_treatment + months(3),
    
    treatment_type == "hormone" &
      is.na(hormone_therapy_end_date1_1)                     ~ NA_POSIXct_,
    treatment_type == "hormone" &
      !is.na(hormone_therapy_end_date1_1)                    ~ hormone_therapy_end_date1_1 + months(3),
    
    str_detect(treatment_type, "chemo") &
      is.na(chemotherapy_end_date1_1)                        ~ NA_POSIXct_,
    str_detect(treatment_type, "chemo") &
      !is.na(chemotherapy_end_date1_1)                       ~ chemotherapy_end_date1_1 + months(3),
    
  )) %>% 
  mutate(reason_for_no_persistent_prolonged_cytopenia_date = case_when(
    treatment_type == "radiation"                            ~ "use radiation start date",
    treatment_type == "hormone" &
      is.na(hormone_therapy_end_date1_1)                     ~ "no hormone end date",
    str_detect(treatment_type, "chemo") &
      is.na(chemotherapy_end_date1_1)                        ~ "no chemo end date"
  )) 

hgb <- HGB_PLT %>%
  filter(lab_nm == "Hemoglobin") %>% 
  right_join(., temp, by = "mrn")  %>% 
  mutate(persistent_low_hgb = case_when(
    lab_date >= (max_persistent_cytopenia_date - months(1)) &
      lab_date <= max_persistent_cytopenia_date &
      lab_result >= 12                                         ~ "No" # 12 is for female
  )) %>% 
  mutate(prolonged_low_hgb = case_when(
    lab_date >= (max_prolonged_cytopenia_date - months(3)) &
      lab_date <= max_prolonged_cytopenia_date &
      lab_result >= 12                                         ~ "No" # 12 is for female
  )) %>% 
  group_by(mrn) %>% 
  fill(persistent_low_hgb, prolonged_low_hgb, .direction = "updown") %>% 
  ungroup() %>% 
  select(mrn, persistent_low_hgb, prolonged_low_hgb, lab_result) %>% 
  distinct(mrn, persistent_low_hgb, prolonged_low_hgb, .keep_all = TRUE) %>% 
  mutate(persistent_low_hgb = case_when(
    persistent_low_hgb == "No"                              ~ "No",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(persistent_low_hgb)                               ~ "Yes"
  )) %>% 
  mutate(prolonged_low_hgb = case_when(
    prolonged_low_hgb == "No"                               ~ "No",
    is.na(lab_result)                                       ~ NA_character_,
    is.na(prolonged_low_hgb)                                ~ "Yes"
  )) %>% 
  select(-lab_result)

plt <- HGB_PLT %>%
  filter(lab_nm == "Platelet Count(k/uL)") %>% 
  right_join(., temp, by = "mrn")  %>% 
  mutate(persistent_low_plt = case_when(
    lab_date >= (max_persistent_cytopenia_date - months(1)) &
      lab_date <= max_persistent_cytopenia_date &
      lab_result >= 150                                      ~ "No"
  )) %>% 
  mutate(prolonged_low_plt = case_when(
    lab_date >= (max_prolonged_cytopenia_date - months(3)) &
      lab_date <= max_prolonged_cytopenia_date &
      lab_result >= 150                                      ~ "No" # 12 is for female
  )) %>% 
  group_by(mrn) %>% 
  fill(persistent_low_plt, prolonged_low_plt, .direction = "updown") %>% 
  ungroup() %>% 
  select(mrn, persistent_low_plt, prolonged_low_plt, lab_result) %>% 
  distinct(mrn, persistent_low_plt, prolonged_low_plt, .keep_all = TRUE) %>% 
  mutate(persistent_low_plt = case_when(
    persistent_low_plt == "No"                              ~ "No",
    is.na(lab_result)                            ~ NA_character_,
    is.na(persistent_low_plt)                               ~ "Yes"
  )) %>% 
  mutate(prolonged_low_plt = case_when(
    prolonged_low_plt == "No"                               ~ "No",
    is.na(lab_result)                            ~ NA_character_,
    is.na(prolonged_low_plt)                                ~ "Yes"
  )) %>% 
  select(-lab_result)

anc <- neutrophil %>%
  right_join(., temp, by = "mrn")  %>% 
  mutate(persistent_low_anc = case_when(
    lab_neutrophil_date >= (max_persistent_cytopenia_date - months(1)) &
      lab_neutrophil_date <= max_persistent_cytopenia_date &
      lab_result_neutrophil >= 1.8                           ~ "No"
  )) %>% 
  mutate(prolonged_low_anc = case_when(
    lab_neutrophil_date >= (max_prolonged_cytopenia_date - months(3)) &
      lab_neutrophil_date <= max_prolonged_cytopenia_date &
      lab_result_neutrophil >= 1.8                           ~ "No"
  )) %>% 
  group_by(mrn) %>% 
  fill(persistent_low_anc, prolonged_low_anc, .direction = "updown") %>% 
  ungroup() %>% 
  select(mrn, persistent_low_anc, prolonged_low_anc, lab_result_neutrophil) %>% 
  distinct(mrn, persistent_low_anc, prolonged_low_anc, .keep_all = TRUE) %>% 
  mutate(persistent_low_anc = case_when(
    persistent_low_anc == "No"                              ~ "No",
    is.na(lab_result_neutrophil)                            ~ NA_character_,
    is.na(persistent_low_anc)                               ~ "Yes"
  )) %>% 
  mutate(prolonged_low_anc = case_when(
    prolonged_low_anc == "No"                               ~ "No",
    is.na(lab_result_neutrophil)                            ~ NA_character_,
    is.na(prolonged_low_anc)                                ~ "Yes"
  )) %>% 
  select(-lab_result_neutrophil)
  
sequenced_patient_data <- sequenced_patient_data %>% 
  left_join(., hgb, by = "mrn") %>% 
  left_join(., plt, by = "mrn") %>% 
  left_join(., anc, by = "mrn") %>% 
  mutate(persistent_cytopenia = case_when(
    persistent_low_hgb == "No" |
      persistent_low_plt == "No" |
      persistent_low_anc == "No"                             ~ "No",
    persistent_low_hgb == "Yes" &
      persistent_low_plt == "Yes" &
      persistent_low_anc == "Yes"                            ~ "Yes",
  )) %>% 
  mutate(prolonged_cytopenia = case_when(
    prolonged_low_hgb == "No" |
      prolonged_low_plt == "No" |
      prolonged_low_anc == "No"                              ~ "No",
    prolonged_low_hgb == "Yes" &
      prolonged_low_plt == "Yes" &
      prolonged_low_anc == "Yes"                             ~ "Yes",
  ))
  
write_csv(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia_",
                                         today(), ".csv"))
write_csv(sequenced_patient_data, paste0(path, "/processed data",
                                         "/Identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia_",
                                         today(), ".csv"))


deids_data <- sequenced_patient_data %>% 
  select(-c(mrn, sample_id, sample_family_id, 
            submitted_ID_for_added_samples, received_ID_for_added_samples,
            Sample_Name_Fastq_file,
            contains("date")
  ))

write_csv(deids_data, 
          paste0(path, 
                 "/processed data", 
                 "/De-identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia_", 
                 today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data", 
                 "/De-identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia_", 
                 today(), ".csv"))


# END NADIR and cytopenia variales----
