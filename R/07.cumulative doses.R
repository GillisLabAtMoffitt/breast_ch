# Import library
library(tidyverse)
library(lubridate)

############################################################ I ### Load data----
path_raw <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "sequential_samples_hossein")
drug_dose <- 
  readxl::read_xlsx(paste0(
    path_raw, 
    "/raw_data/Chemo_patient_MRN_chemotherapy_drug_04.02.2025 1.xlsx"),
    sheet = "Normalized") %>% 
  janitor::clean_names()

sample_dates <- 
  readxl::read_xlsx(paste0(
    path_raw, 
    "/raw_data/Chemo_patient_MRN_chemotherapy_drug_04.02.2025 1.xlsx"),
    sheet = "Chemo_patient_MRN_chemotherapy_") %>% 
  janitor::clean_names()

drug_class <- 
  readxl::read_xlsx(paste0(
    path_raw, 
    "/raw_data/Bolton_Supp Tables.xlsx"),
    sheet = "Sup. Table 1", skip = 1) %>% 
  janitor::clean_names()

path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- 
  read_csv(paste0(# path_save, 
    here::here(), 
    "/processed data",
    "/Identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia-NADIR_2025-03-31.csv"))


############################################################ II ### Cumulative doses----
sample_dates <- sample_dates %>% 
  select(mrn, pre_sample_collection_date, last_post_sample_collection_date) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  mutate_at(c("pre_sample_collection_date", "last_post_sample_collection_date"), ~ as.Date(., origin = "1899-12-30")) %>% 
  left_join(., sequenced_patient_data %>% 
              # select(mrn, treatment_type, seqsample_date_chemo, seqsample_date_rad, 
              #        seqsample_date_Achemo, seqsample_date_hormone) %>% 
              group_by(mrn) %>% 
              fill(seqsample_date_chemo, seqsample_date_rad, 
                   seqsample_date_Achemo, seqsample_date_hormone, .direction = "updown") %>% 
              distinct(mrn, .keep_all = TRUE) %>% 
              mutate(date_of_first_seq_sample = case_when(
                treatment_type == "chemo"              ~ seqsample_date_chemo,
                treatment_type == "chemorad"           ~ seqsample_date_Achemo,
                treatment_type == "hormone"            ~ seqsample_date_hormone,
                treatment_type == "radiation"          ~ seqsample_date_rad
              )) %>% 
              select(mrn, date_of_first_seq_sample),
            by = "mrn") %>% 
  mutate(mrn = as.character(mrn))
  
drug_dose1 <- drug_dose %>% 
  mutate(mrn = as.character(mrn)) %>% 
  left_join(., sample_dates, by = "mrn") %>% 
  mutate(pre_to_first_sample = case_when(
    date >= pre_sample_collection_date &
      date <= date_of_first_seq_sample                  ~ "Yes"
  )) %>% 
  mutate(pre_to_last_sample = case_when(
    date >= pre_sample_collection_date &
      date <= last_post_sample_collection_date          ~ "Yes"
  )) %>% 
  # fix drug name and add class
  mutate(chemotherapy_drug = case_when(
    chemotherapy_drug == "5 fu"                         ~ "fluorouracil",
    chemotherapy_drug == "adriamycin"                   ~ "doxorubicin",
    chemotherapy_drug == "abraxane (form of taxol)"     ~ "abraxane",
    TRUE                                                ~ chemotherapy_drug
  )) %>% 
  # add 
  left_join(., drug_class, 
            by = c("chemotherapy_drug" = "drug_name")) %>% 
  # calculate cumulative dose per patient per drug
  group_by(mrn, chemotherapy_drug, pre_to_first_sample) %>% 
  mutate(cumulative_dose_pre_firstseqsample = case_when(
    pre_to_first_sample == "Yes"                        ~ sum(dose_mg_per_m2)
  )) %>% 
  group_by(mrn, chemotherapy_drug, pre_to_last_sample) %>% 
  mutate(cumulative_dose_pre_lastseqsample = case_when(
    pre_to_last_sample == "Yes"                         ~ sum(dose_mg_per_m2)
  )) %>% 
  ungroup() %>% 
  select(mrn, chemotherapy_drug, narrow_drug_class_cytotoxic_only, 
         cumulative_dose_pre_firstseqsample, cumulative_dose_pre_lastseqsample) %>% 
  arrange(mrn, chemotherapy_drug, 
          cumulative_dose_pre_firstseqsample, cumulative_dose_pre_lastseqsample) %>% 
  distinct(mrn, chemotherapy_drug, .keep_all = TRUE) %>% 
  # create tertile categories
  group_by(chemotherapy_drug) %>% 
  # mutate(n = n()) %>%
  # mutate(tertile_pre_firstseqsample = ntile(cumulative_dose_pre_firstseqsample, 3))
  mutate(tertile_pre_firstseqsample = case_when(
    n() > 2                                             ~ ntile(cumulative_dose_pre_firstseqsample, 3),
    n() == 1 &
      chemotherapy_drug == "abraxane"                   ~ 1,
    n() == 2 &
      chemotherapy_drug %in% 
      c("epirubicin", "fluorouracil")                   ~ 2
  )) %>% 
  mutate(tertile_pre_lastseqsample = case_when(
    n() > 2                                             ~ ntile(cumulative_dose_pre_lastseqsample, 3),
    n() == 1 &
      chemotherapy_drug == "abraxane"                   ~ 1,
    n() == 2 &
      chemotherapy_drug %in% 
      c("epirubicin", "fluorouracil")                   ~ 2
  )) %>% 
  # mutate(tertile_cum_drug_cat = case_when(
  #   n() > 1 &
  #     ntile(cumulative_dose, 3) == 1                    ~ "Low",
  #   n() > 1 &
  #     ntile(cumulative_dose, 3) == 2                    ~ "Medium",
  #   n() > 1 &
  #     ntile(cumulative_dose, 3) == 3                    ~ "High",
  # )) %>% 
  
  # calculate score - sum scores for each drug in a specific drug class
  group_by(mrn, narrow_drug_class_cytotoxic_only) %>% 
  # mutate(n = n()) %>%
  mutate(drug_class_score_pre_firstseqsample = case_when(
    !is.na(tertile_pre_firstseqsample)                  ~ sum(tertile_pre_firstseqsample, na.rm = TRUE)
  )) %>% 
  mutate(drug_class_score_pre_lastseqsample = case_when(
    !is.na(tertile_pre_lastseqsample)                   ~ sum(tertile_pre_lastseqsample, na.rm = TRUE)
  )) %>% 
  distinct(mrn, narrow_drug_class_cytotoxic_only, .keep_all = TRUE) %>% 
  select(-chemotherapy_drug) %>% 
  # group_by(narrow_drug_class_cytotoxic_only) %>%
  # mutate(n = n()) %>% 
  # Bolton summed score
  group_by(mrn, narrow_drug_class_cytotoxic_only) %>% 
  mutate(sum_bolton_drug_score_pre_firstseqsample_per_class_patient = case_when(
    !is.na(drug_class_score_pre_firstseqsample)         ~ sum(drug_class_score_pre_firstseqsample, na.rm = TRUE)
  )) %>% 
  mutate(sum_bolton_drug_score_pre_lastseqsample_per_class_patient = case_when(
    !is.na(drug_class_score_pre_lastseqsample)          ~ sum(drug_class_score_pre_lastseqsample, na.rm = TRUE)
  )) %>% 
  # Patient's
  group_by(mrn) %>% 
  mutate(sum_all_drug_score_pre_firstseqsample_per_patient = case_when(
    !is.na(drug_class_score_pre_firstseqsample)         ~ sum(drug_class_score_pre_firstseqsample, na.rm = TRUE)
  )) %>% 
  mutate(sum_all_drug_score_pre_lastseqsample_per_patient = case_when(
    !is.na(drug_class_score_pre_lastseqsample)          ~ sum(drug_class_score_pre_lastseqsample, na.rm = TRUE)
  )) %>% 
  ungroup() %>% 
  # # use cut() instead of ntile() to attribute the ties in values
  mutate(tertile_patient_pre_firstseqsample = cut(sum_all_drug_score_pre_firstseqsample_per_patient,
                             breaks = quantile(sum_all_drug_score_pre_firstseqsample_per_patient, 
                                               probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                             include.lowest = TRUE,
                             labels = c("Low", "Medium", "High"))) %>%
  mutate(tertile_patient_pre_lastseqsample = cut(sum_all_drug_score_pre_lastseqsample_per_patient,
                                                breaks = quantile(sum_all_drug_score_pre_lastseqsample_per_patient, 
                                                                  probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                                include.lowest = TRUE,
                                                labels = c("Low", "Medium", "High"))) %>% 
  full_join(sequenced_patient_data %>% 
              select(mrn, deidentified_patient_id) %>% 
              distinct() %>% 
              mutate(mrn = as.character(mrn)), 
            ., 
            by = "mrn") %>% 
  select(mrn, deidentified_patient_id, narrow_drug_class_cytotoxic_only, 
         sum_bolton_drug_score_pre_firstseqsample_per_class_patient, 
         sum_bolton_drug_score_pre_lastseqsample_per_class_patient,
         sum_all_drug_score_pre_firstseqsample_per_patient,
         tertile_patient_pre_firstseqsample,
         sum_all_drug_score_pre_lastseqsample_per_patient,
         tertile_patient_pre_lastseqsample
         ) %>% 
  filter(!is.na(narrow_drug_class_cytotoxic_only))
  

write_csv(drug_dose1, 
          paste0(here::here(), 
                 "/processed data",
                 "/Cumulative dose score_",
                 today(), ".csv"))
write_rds(drug_dose1, 
          paste0(here::here(), 
                 "/processed data",
                 "/Cumulative dose score_",
                 today(), ".rds"))
write_csv(drug_dose1, 
          paste0(path_save, 
                 "/processed data",
                 "/Cumulative dose score_",
                 today(), ".csv"))
write_csv(drug_dose1 %>% 
            select(-mrn), 
          paste0(path_save, 
                 "/processed data",
                 "/De-identified Cumulative dose score_",
                 today(), ".csv"))










