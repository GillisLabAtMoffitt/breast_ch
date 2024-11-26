# Code "add-post" samples

# Import library
library(tidyverse)

# Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                 "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read.csv(paste0(#path_save,
  here::here(),
  "/processed data/Identified breast data with sequenced sequential sample and clinical_2024-09-30.csv")) %>% 
  mutate(mrn = as.character(mrn))
sequenced_patient_data <- read_rds(paste0(#path_save,
  here::here(),
  "/processed data/Identified breast data with sequenced sequential sample and clinical_2024-09-30.rds"))

path_raw <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "sequential_samples_hossein")
chart_dat <- 
  readxl::read_xlsx(paste0(#path_raw, 
    here::here(),
    "/chart_reviewed/Breast_chart reviews template_11252024.xlsx"),
                    sheet = "Treatment_CBCs", na = c("NA", "UNK", "ND")) %>% 
  janitor::clean_names()


# Cleaning----
chart_dat <- chart_dat %>% 
  # Remove duplicated rows which appears between my output file and what Danny received
  # filter(is.na(duplicated_rows)) %>% # this variable is wrong
  group_by(mrn, sample_name_secondary, collection_date) %>% 
  fill(everything(), .direction = "updown") %>% 
  ungroup() %>% 
  distinct(mrn, sample_name_secondary, collection_date, .keep_all = TRUE) %>% 
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
  # mutate(#progression_comment = date_of_progression, .after = date_of_progression,
  #        date_of_progression = as.Date(as.numeric(date_of_progression),
  #                                      # format = "%Y-%m-%d"), 
  #                                      origin = "1899-12-30"
  #        # ),
  #        # progression_comment = case_when(
  #        #   is.na(date_of_progression)               ~ progression_comment # No need anymore
  #        )) %>% 
  # mutate(#last_followup_comment = date_of_last_followup, .after = date_of_last_followup,
  #        date_of_last_followup = as.Date(as.numeric(date_of_last_followup),
  #                                        # format = "%Y-%m-%d"), 
  #                                        origin = "1899-12-30"
  #        # ),
  #        # last_followup_comment = case_when(
  #        #   is.na(date_of_last_followup)             ~ last_followup_comment # No need anymore
  #        )) %>% 
  select(-c(x0, x63))
  

chart_dat1 <- chart_dat %>% 
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
  select(-boost_dose_c_gy_1) %>% 
  full_join(.,
            chart_dat1 %>% select(-collection_date), 
            by = c("mrn", 
                   "Sample.Name..Secondary." = "sample_name_secondary")) %>% 
  full_join(.,
            chart_dat %>% 
              select(mrn, type_of_radiation : boost_dose_c_gy_1) %>% 
              distinct(), 
            by = c("mrn")) %>% 
  mutate(type_of_radiation_1 = type_of_radiation, .after = radiation_start_date_1) %>% 
  mutate(total_dose_1 = total_dose, .after = type_of_radiation_1) %>% 
  mutate(boost_dose_c_gy = boost_dose_c_gy_1, .after = total_dose_1) %>% 
  select(mrn : specimen_collection_date, 
         age_at_sample : treatment_type, sequence_sample_vs_treatment, treatment_received,
         everything(), -type_of_radiation, -total_dose) %>% 
  mutate(discrepancies_old_new_sample_sequence = case_when(
    treatment_type == "chemo" &
      str_detect(sequence_sample_vs_treatment, "after_hormone|after_radiation")      ~ "need to be reclassified/remove",
    treatment_type == "chemorad" &
      str_detect(sequence_sample_vs_treatment, "after_hormone")                      ~ "need to be reclassified/remove",
    
  ), .after = treatment_type)



# # Save----
# write_csv(sequenced_patient_data, 
#           paste0(path_save, 
#                  "/processed data/Identified breast data with re-classified sequenced sequential sample and clinical_", 
#                  today(), ".csv"))
# write_csv(sequenced_patient_data, 
#           paste0("processed data/Identified breast data with re-classified sequenced sequential sample and clinical_", 
#                  today(), ".csv"))
# write_rds(sequenced_patient_data, 
#           paste0("processed data/Identified breast data with re-classified sequenced sequential sample and clinical_", 
#                  today(), ".rds"))
# 
# deids_data <- sequenced_patient_data %>% 
#   select(-c(mrn, sample_id, sample_family_id, Sample.Name..Secondary.,
#             submitted_ID_for_added_samples, received_ID_for_added_samples,
#             Sample_Name_Fastq_file,
#             contains("date")
#   ))
# 
# write_csv(deids_data, 
#           paste0(path_save, 
#                  "/processed data/De-identified breast data with re-classified sequenced sequential sample_", today(), ".csv"))
# write_csv(deids_data, 
#           paste0("processed data/De-identified breast data with re-classified sequenced sequential sample_", 
#                  today(), ".csv"))

# END re-classified samples ----

# Update clinical with chart review----
sequenced_patient_data <- chart_dat %>% 
  # Select 1 clinical row for each patient
  select(mrn, first_treatment_h_end_date,
         smoking_status : oncotype_dx,
         dead_or_alive,
         date_of_last_followup, 
         progression_type, date_of_progression) %>% 
  group_by(mrn) %>%
  fill(everything(), .direction = "updown") %>%
  ungroup() %>%
  distinct(mrn, .keep_all = TRUE) %>% 
  mutate_at(c("dead_or_alive"), ~ str_to_sentence(.)) %>% 
  left_join(sequenced_patient_data %>% 
              # Remove initially created variables before to add updated chart reviewed var
              select(-c(
                smoking_status, new_vital_status,
                contains("progesterone"), contains("estrogen"),
                starts_with("os_"), type_of_first_recurrence_desc,
                ends_with("_results"), starts_with("her2"), 
                ER_PR_status, ER_PR_HER_status, 
                contains("oncotype_dx")
              )), 
            ., 
            by = "mrn") %>% 
  # Update OS and PFS
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
  
  mutate(pfs_event = case_when(
    is.na(progression_type)             ~ 0,
    !is.na(progression_type)            ~ 1
  ), .after = progression_type) %>% 
  mutate(pfs_date = coalesce(date_of_progression, date_of_last_followup), 
         .after = pfs_event) %>% 
  mutate(pfs_time_months = interval(start = date_of_first_corresponding_treatment, end = pfs_date)/
           duration(n = 1, unit = "months"), 
         .after = pfs_event) %>% 
  # update chemo length as Nancy changed a end date
  mutate(chemo_length_indays = case_when(
    treatment_type == "chemo" |
      treatment_type == "chemorad"         ~ interval(start = chemotherapy_start_date_1,
                                                      end = first_treatment_h_end_date)/
      duration(n = 1, units = "days")
  )) %>% 
  mutate(chart_reviewed_chemotherapy_end_date = first_treatment_h_end_date,
         .after = chemotherapy_end_date1_1) %>% 
  select(-c(first_treatment_h_end_date, chemotherapy_end_date1_1))

write_csv(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data with re-classified sequenced sequential sample and chart review clinical_",
                                         today(), ".csv"))
write_rds(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data with re-classified sequenced sequential sample and chart review clinical_",
                                         today(), ".rds"))
write_csv(sequenced_patient_data, paste0(path_save, "/processed data",
                                         "/Identified breast data with re-classified sequenced sequential sample and chart review clinical_",
                                         today(), ".csv"))

# Save de-identified data
deids_data <- sequenced_patient_data %>% 
  select(-c(mrn, sample_id, sample_family_id, Sample.Name..Secondary.,
            submitted_ID_for_added_samples, received_ID_for_added_samples,
            Sample_Name_Fastq_file,
            contains("date")
  ))

write_csv(deids_data, 
          paste0(path_save, 
                 "/processed data", 
                 "/De-identified breast data with re-classified sequenced sequential sample and chart review clinical_", 
                 today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data", 
                 "/De-identified breast data with re-classified sequenced sequential sample and chart review clinical_", 
                 today(), ".csv"))


# Patient without pre-sample CBC
a <- sequenced_patient_data %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, overall_neutroplil_lab_result_at_presample, starts_with("lab_result_") & ends_with("_at_presample")) %>%
  filter(if_all(any_of(c(contains("lab_result_"))), ~ is.na(.)))

c <- sequenced_patient_data %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, 
         lab_result_wbc_at_presample, lab_result_Hemoglobin_at_presample,
         lab_result_Platelet.Count.k.uL._at_presample,
         lab_result_mcv_at_presample, lab_result_rdw_at_presample,
         overall_neutroplil_lab_result_at_presample
         ) %>%
  filter(if_any(any_of(c(contains("lab_result_"))), ~ is.na(.))) %>% 
  filter(!str_detect(mrn, paste0(a$mrn, collapse = "|")
                     )) %>% 
  arrange(overall_neutroplil_lab_result_at_presample)


b <- sequenced_patient_data %>% select(starts_with("lab_result_") & ends_with("_at_presample"))
         
write_csv(a, 
          "breast patients with no pre sample cbc at all.csv")
write_csv(c, 
          "breast patients with at least 1 missing pre sample cbc arrange by anc.csv")


# END Update clinical with chart review----

