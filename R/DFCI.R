library(tidyverse)

############################################################ I ### Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read_rds(paste0(
  here::here(),
  "/processed data",
  "/Identified breast data_2025-06-30.rds"))

path_raw <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "Dana Farber breast data")
patient_data <- 
  read_csv(
    paste0(path_raw, "/raw data",
           "/DFCI.CH.patients2025-07-08.csv")) %>% 
  janitor::clean_names()

chemo_data <- 
  read_csv(
    paste0(path_raw, "/raw data",
           "/DFCI.CH.chemo2025-07-08.csv")) %>% 
  janitor::clean_names()

rad_data <- 
  read_csv(
    paste0(path_raw, "/raw data",
           "/DFCI.CH.radiation2025-07-08.csv")) %>% 
  janitor::clean_names()

path_raw2 <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "sequential_samples_hossein")

drug_class <- 
  readxl::read_xlsx(paste0(
    path_raw2, 
    "/raw_data/Bolton_Supp Tables.xlsx"),
    sheet = "Sup. Table 1", skip = 1) %>% 
  janitor::clean_names()


############################################################ II ### Clean data----
patient_data <- patient_data %>% 
  rename(deidentified_patient_id = id, 
         age_at_sample = age,
         smoking_status = smoking,
         treatment_type = cohort
         # interval_prepost_sample_chemo = time_between_draws
  ) %>% 
  mutate(race = str_to_sentence(race),
         race = case_when(
           race == "Caucasian"                             ~ "White",
           race == "Black or african american"             ~ "Black /african american",
           race == "Asian or pacific islander"             ~ "Asian",
           race == "American indian, aleutian, eskimo"     ~ "American indian or alaska native",
           race == "Other"                                 ~ "Other race",
           # keep their unknown as unknown
           race == "Unknown"                               ~ "Unknown",
         )) %>% 
  group_by(deidentified_patient_id) %>% 
  mutate(tnm_stage = max(clinical_stage, pathological_stage), .after = pathological_stage) %>% 
  ungroup() %>% 
  mutate(tumor_subtype = case_when(
    tumor_subtype == "TNBC"             ~ "HR-/HER2-",
    TRUE                                ~ tumor_subtype
  )) %>% 
  separate_wider_delim(cols = tumor_subtype, delim = "/",
                       names = c("hormone_status", "her2"), 
                       too_few = "align_start", too_many = "error", 
                       cols_remove = FALSE) %>% 
  mutate(er = case_when(
    hormone_status == "HR-"            ~ "neg",
    hormone_status == "HR+"            ~ "pos"
  )) %>% 
  mutate(pr = case_when(
    hormone_status == "HR-"            ~ "neg",
    hormone_status == "HR+"            ~ "pos"
  )) %>% 
  mutate(her2 = case_when(
    her2 == "HER2-"                    ~ "neg",
    her2 == "HER2+"                    ~ "pos"
  ))



chemo_data1 <- chemo_data %>% 
  rename(deidentified_patient_id = id) %>% 
  group_by(deidentified_patient_id, regimen_name) %>%
  summarise_at(vars(drug_name, number_of_cycles), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  group_by(deidentified_patient_id) %>%
  summarise_at(vars(regimen_name, drug_name, number_of_cycles), 
               str_c, collapse = " AND ") %>% 
  ungroup() %>% 
  select(deidentified_patient_id, chemotherapy_drug_1 = drug_name)

rad_data1 <- rad_data %>% 
  select(id, had_rad = rt_between_t1_t2)



patient_data1 <- patient_data %>% 
  # Update the 1 patient who only received targeted therapy
  left_join(., chemo_data1, 
            by = "deidentified_patient_id") %>% 
  mutate_at(c("treatment_type"), ~ case_when(
    chemotherapy_drug_1 == "T-DM1"                  ~ "B - no chemo",
    TRUE                                            ~ as.character(.)
  )) %>% 
  # Add radiation info
  left_join(., rad_data1, 
            by = c("deidentified_patient_id" = "id")) %>% 
  
  mutate(treatment_type = case_when(
    had_rad == "Yes" &
      treatment_type == "A - chemo"             ~ "chemorad",
    treatment_type == "A - chemo"               ~ "chemo",
    had_rad == "Yes" &
      treatment_type == "B - no chemo"          ~ "radiation",
    treatment_type == "B - no chemo"            ~ "hormone"#,
    # had_rad == "Yes"                            ~ "radiation"
  )) %>% 
  mutate(treatment_received = treatment_type) %>% 

  select(deidentified_patient_id, age_at_sample, race,
         smoking_status, prior_cancer, prior_radiation,
         er, pr, her2, 
         histology, tnm_stage, had_rad,
         treatment_type, treatment_received, time_between_draws) %>% 
  
  mutate(interval_prepost_sample_rad = case_when(
    treatment_type == "radiation"      ~ time_between_draws
  )) %>% 
  mutate(interval_prepost_sample_chemo = case_when(
    treatment_type == "chemo"          ~ time_between_draws
  )) %>% 
  mutate(had_chemo = case_when(
    treatment_type == "chemo"          ~ "Yes",
    treatment_type == "chemorad"       ~ "Yes"
  )) %>% 
  mutate(interval_prepost_sample_hormone = case_when(
    treatment_type == "hormone"        ~ time_between_draws
  )) %>% 
  mutate(had_hormone = case_when(
    treatment_type == "hormone"        ~ "Yes"
  )) %>% 
  
  mutate(interval_prepost_sample_chemorad = case_when(
    treatment_type == "chemorad"       ~ time_between_draws
  )) %>% 
  mutate(interval_prepost_sample_chemo_rad = interval_prepost_sample_chemorad) %>% 
  mutate(had_chemo_rad = case_when(
    treatment_type == "chemorad"       ~ "Yes"
  )) %>% 
  select(-time_between_draws)


############################################################ III ### Merge with data format----
patient_data1 <- sequenced_patient_data %>% 
  slice(1) %>% 
  bind_rows(., patient_data1) %>% 
  slice(-1)

write_csv(patient_data1, paste0(here::here(), "/processed data",
                                         "/DFCI clinical breast data_",
                                         today(), ".csv"))
write_csv(patient_data1, paste0(path_save, "/processed data",
                                         "/DFCI clinical breast data_",
                                         today(), ".csv"))
write_rds(patient_data1, paste0(path_save, "/processed data",
                                         "/DFCI clinical breast data_",
                                         today(), ".rds"))
write_rds(patient_data1, paste0(here::here(), "/processed data",
                                         "/DFCI clinical breast data_",
                                         today(), ".rds"))


deids_data <- patient_data1 %>% 
  select(-c(mrn, sample_id, sample_family_id, sample_name_secondary, party_id,
            submitted_ID_for_added_samples, received_ID_for_added_samples,
            Sample_Name_Fastq_file,
            contains("date"),
            blood_bf_chemo_rad : sequenced_samples
  ))

write_csv(deids_data, 
          paste0(path_save, 
                 "/processed data", 
                 "/Limited DFCI clinical breast data_", 
                 today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data", 
                 "/Limited DFCI clinical breast data_", 
                 today(), ".csv"))


############################################################ IV ### Cumulative doses----
chemo_data2 <- chemo_data %>% 
  select(deidentified_patient_id = id,
         # chemotherapy_drug_1 = 
         drug_name,
         number_of_cycles) %>% 
  left_join(., drug_class, 
            by = "drug_name") %>% 
  filter(!is.na(narrow_drug_class_cytotoxic_only)) %>% 
  
  
  # calculate cumulative dose per patient per drug
  group_by(deidentified_patient_id, drug_name) %>% 
  mutate(cumulative_dose_pre_firstseqsample = sum(as.numeric(number_of_cycles))
  ) %>% 
  ungroup() %>% 
  select(deidentified_patient_id, drug_name, narrow_drug_class_cytotoxic_only, 
         cumulative_dose_pre_firstseqsample) %>% 
  arrange(deidentified_patient_id, drug_name, 
          cumulative_dose_pre_firstseqsample) %>% 
  distinct(deidentified_patient_id, drug_name, .keep_all = TRUE) %>% 
  # create tertile categories
  group_by(drug_name) %>% 
  # mutate(n = n()) %>% 
  # mutate(tertile_pre_firstseqsample = ntile(cumulative_dose_pre_firstseqsample, 3))
  mutate(tertile_pre_firstseqsample = case_when(
    n() > 2                                             ~ ntile(cumulative_dose_pre_firstseqsample, 3),
    n() == 1 &
      drug_name == "abraxane"                           ~ 1,
    n() == 2 &
      drug_name %in% 
      c("methotrexate", "fluorouracil")                 ~ 2
  )) %>% 
  group_by(deidentified_patient_id, narrow_drug_class_cytotoxic_only) %>% 
  # mutate(n = n()) %>%
  mutate(drug_class_score_pre_firstseqsample = case_when(
    !is.na(tertile_pre_firstseqsample)                  ~ sum(tertile_pre_firstseqsample, na.rm = TRUE)
  )) %>% 

  distinct(deidentified_patient_id, narrow_drug_class_cytotoxic_only, .keep_all = TRUE) %>% 
  select(-drug_name) %>% 
  group_by(deidentified_patient_id, narrow_drug_class_cytotoxic_only) %>% 
  mutate(sum_bolton_drug_score_pre_firstseqsample_per_class_patient = case_when(
    !is.na(drug_class_score_pre_firstseqsample)         ~ sum(drug_class_score_pre_firstseqsample, na.rm = TRUE)
  )) %>% 

  # Patient's
  group_by(deidentified_patient_id) %>% 
  mutate(sum_all_drug_score_pre_firstseqsample_per_patient = case_when(
    !is.na(drug_class_score_pre_firstseqsample)         ~ sum(drug_class_score_pre_firstseqsample, na.rm = TRUE)
  )) %>% 
  ungroup() %>% 
  # # use cut() instead of ntile() to attribute the ties in values
  mutate(tertile_patient_pre_firstseqsample = cut(sum_all_drug_score_pre_firstseqsample_per_patient,
                                                  breaks = quantile(sum_all_drug_score_pre_firstseqsample_per_patient, 
                                                                    probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                                  include.lowest = TRUE,
                                                  labels = c("Low", "Medium", "High"))) %>% 
  select(deidentified_patient_id, narrow_drug_class_cytotoxic_only, 
         sum_bolton_drug_score_pre_firstseqsample_per_class_patient, 
         sum_all_drug_score_pre_firstseqsample_per_patient,
         tertile_patient_pre_firstseqsample
  ) %>% 
  filter(!is.na(narrow_drug_class_cytotoxic_only))

write_csv(chemo_data2, 
          paste0(here::here(), 
                 "/processed data",
                 "/DFCI Cumulative dose score_",
                 today(), ".csv"))
write_rds(chemo_data2, 
          paste0(here::here(), 
                 "/processed data",
                 "/DFCI Cumulative dose score_",
                 today(), ".rds"))
write_csv(chemo_data2, 
          paste0(path_save, 
                 "/processed data",
                 "/DFCI Cumulative dose score_",
                 today(), ".csv"))


