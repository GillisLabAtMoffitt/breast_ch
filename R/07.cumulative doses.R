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
  select(mrn, chemotherapy_drug, general_drug_class,
         cumulative_dose_pre_firstseqsample, cumulative_dose_pre_lastseqsample) %>% 
  arrange(mrn, chemotherapy_drug, 
          cumulative_dose_pre_firstseqsample, cumulative_dose_pre_lastseqsample) %>% 
  distinct(mrn, chemotherapy_drug, .keep_all = TRUE) %>% 
  # create tertile categories
  group_by(chemotherapy_drug) %>% 
  # mutate(n = n()) %>% 
  mutate(tertile_pre_firstseqsample = case_when(
    n() > 1                                             ~ ntile(cumulative_dose_pre_firstseqsample, 3)
  )) %>% 
  mutate(tertile_pre_lastseqsample = case_when(
    n() > 1                                             ~ ntile(cumulative_dose_pre_lastseqsample, 3)
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
  group_by(mrn, general_drug_class) %>% 
  mutate(drug_class_score_pre_firstseqsample = sum(tertile_pre_firstseqsample, na.rm = TRUE)) %>% 
  mutate(drug_class_score_pre_lastseqsample = sum(tertile_pre_lastseqsample, na.rm = TRUE)) %>% 
  distinct(mrn, general_drug_class, .keep_all = TRUE) %>% 
  # group_by(general_drug_class) %>% 
  # mutate(n = n()) %>% 
  # mutate(tertile_class = case_when(
  #   n() > 1                                             ~ ntile(drug_class_score, 3)
  # )) %>%
  # # calculate score - tertile for a general_drug_class
  # mutate(drug_class_score_cat = case_when(
  #   n() > 1 &
  #     ntile(floor(drug_class_score), 3) == 1                   ~ "Low",
  #   n() > 1 &
  #     ntile(floor(drug_class_score), 3) == 2                   ~ "Medium",
  #   n() > 1 &
  #     ntile(floor(drug_class_score), 3) == 3                   ~ "High",
  # )) %>% 
  ungroup() %>% 
  # use cut() instead of ntile() to get the ties in values
  mutate(tertile_class_pre_firstseqsample = cut(drug_class_score_pre_firstseqsample,
                             breaks = quantile(drug_class_score_pre_firstseqsample, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                             include.lowest = TRUE,
                             labels = c("Low", "Medium", "High"))) %>% 
  mutate(tertile_class_pre_lastseqsample = cut(drug_class_score_pre_lastseqsample,
                                                breaks = quantile(drug_class_score_pre_lastseqsample, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                                include.lowest = TRUE,
                                                labels = c("Low", "Medium", "High"))) %>% 
  select(mrn, general_drug_class, tertile_class_pre_firstseqsample, 
         tertile_class_pre_lastseqsample, drug_class_score_pre_firstseqsample, drug_class_score_pre_lastseqsample)


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










