# Code "add-post" samples

# Import library
library(tidyverse)

# Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                 "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read.csv(paste0(#path_save,
  here::here(),
  "/processed data/Identified breast data with sequenced sequential sample and clinical_2024-09-30.csv"))
sequenced_patient_data <- read_rds(paste0(#path_save,
  here::here(),
  "/processed data/Identified breast data with sequenced sequential sample and clinical_2024-09-30.rds"))

path_raw <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "sequential_samples_hossein")
chart_dat <- 
  readxl::read_xlsx(paste0(#path_raw, 
    here::here(),
    "/chart_reviewed/Breast_chart reviews template_11052024.xlsx"),
                    sheet = "Treatment_CBCs", na = c("NA", "UNK")) %>% 
  janitor::clean_names()


# Cleaning----
chart_dat <- chart_dat %>% 
  filter(is.na(duplicated_rows)) %>% 
  mutate(first_treatment_h_end_date = as.Date(as.numeric(first_treatment_h_end_date),
                                              # format = "%Y-%m-%d"), 
                                              origin = "1899-12-30"
                                              )) %>% 
  mutate(hormone_end_date_current = case_when(
    hormone_end_date1_1 == "Current"                ~ "Currently use"
  ), .after = hormone_end_date1_1) %>% 
  mutate(hormone_end_date1_1 = as.Date(as.numeric(hormone_end_date1_1),
                                              # format = "%Y-%m-%d"), 
                                              origin = "1899-12-30"
  )) %>% 
  mutate(hormone_end_date1_1 = case_when(
    !is.na(hormone_end_date1_1)                     ~ hormone_end_date1_1,
    hormone_end_date_current == "Currently use"     ~ today() + months(1)
  )) %>% 
  mutate(progression_comment = date_of_progression, .after = date_of_progression,
         date_of_progression = as.Date(as.numeric(date_of_progression),
                                       # format = "%Y-%m-%d"), 
                                       origin = "1899-12-30"
         ),
         progression_comment = case_when(
           is.na(date_of_progression)               ~ progression_comment
         )) %>% 
  mutate(last_followup_comment = date_of_last_followup, .after = date_of_last_followup,
         date_of_last_followup = as.Date(as.numeric(date_of_last_followup),
                                         # format = "%Y-%m-%d"), 
                                         origin = "1899-12-30"
         ),
         last_followup_comment = case_when(
           is.na(date_of_last_followup)             ~ last_followup_comment
         ))
  

chart_dat <- chart_dat %>% 
  pivot_longer(cols = c(chemo1_start_date_1, chemo2_start,
                        radiation_start_date_1, hormone_start_date_1), 
               names_to = "treatment_type", 
               values_to = "treatment_sequence_date", 
               values_drop_na = TRUE) %>% 
  mutate(treatment_type = str_extract(treatment_type, "chemo1|chemo2|radiation|hormone")) %>% 
  arrange(mrn, collection_date, treatment_sequence_date) %>% 
  select(mrn, sample_name_secondary, collection_date, 
         treatment_sequence_date,
         treatment_type) %>% 
  mutate(sample_vs_treatment = case_when(
    collection_date <= treatment_sequence_date          ~ "before",
    collection_date > treatment_sequence_date           ~ "after"
  )) %>% 
  unite(sequence, c(sample_vs_treatment, treatment_type), remove = FALSE) %>% 
  group_by(mrn, sample_name_secondary) %>% 
  mutate(treatment_number = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(mrn, sample_name_secondary, collection_date), 
              names_from = treatment_number, 
              values_from = c(sequence, treatment_type)) %>% 
  unite(treatment_received, starts_with("treatment_type"), na.rm = TRUE, remove = TRUE) %>% 
  unite(sequence_sample_vs_treatment, starts_with("sequence"), sep = ", ", na.rm = TRUE, remove = TRUE) %>% 
  mutate(sequence_sample_vs_treatment = case_when(
    !str_detect(sequence_sample_vs_treatment, "after")             ~ "pre",
    TRUE               ~ sequence_sample_vs_treatment
  ))
  

# Join new sample classification with data----
sequenced_patient_data <- sequenced_patient_data %>% 
  mutate(specimen_collection_date = as.POSIXct(specimen_collection_date, origin = "1899-12-30")) %>% 
  full_join(.,
            chart_dat, 
            by = c("mrn", 
                   "Sample.Name..Secondary." = "sample_name_secondary")) %>% 
  select(mrn : specimen_collection_date, collection_date_danny = collection_date,
         age_at_sample : treatment_type, sequence_sample_vs_treatment, treatment_received,
         everything()) %>% 
  mutate(discrepancies_old_new_sample_sequence = case_when(
    treatment_type == "chemo" &
      str_detect(sequence_sample_vs_treatment, "after_hormone|after_radiation")      ~ "need to be reclassified/remove",
    treatment_type == "chemorad" &
      str_detect(sequence_sample_vs_treatment, "after_hormone")                      ~ "need to be reclassified/remove",
    
  ), .after = treatment_type)



# Save----
write_csv(sequenced_patient_data, 
          paste0(path_save, 
                 "/processed data/Identified breast data with re-classified sequenced sequential sample and clinical_", 
                 today(), ".csv"))
write_csv(sequenced_patient_data, 
          paste0("processed data/Identified breast data with re-classified sequenced sequential sample and clinical_", 
                 today(), ".csv"))
write_rds(sequenced_patient_data, 
          paste0("processed data/Identified breast data with re-classified sequenced sequential sample and clinical_", 
                 today(), ".rds"))

deids_data <- sequenced_patient_data %>% 
  select(-c(mrn, sample_id, sample_family_id, Sample.Name..Secondary.,
            submitted_ID_for_added_samples, received_ID_for_added_samples,
            Sample_Name_Fastq_file,
            contains("date")
  ))

write_csv(deids_data, 
          paste0(path_save, 
                 "/processed data/De-identified breast data with re-classified sequenced sequential sample_", today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data/De-identified breast data with re-classified sequenced sequential sample_", 
                 today(), ".csv"))

# END re-classified samples ----



