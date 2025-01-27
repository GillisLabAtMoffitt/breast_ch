# Import library
library(tidyverse)
library(lubridate)

############################################################ I ### Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

blood_patients <- 
  read_rds(paste0(here::here(), "/blood_patients2.rds"))

later_blood_patients <-
  read_rds(paste0(here::here(), "/blood_patients_06272023.rds"))

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
  read_csv(paste0(here::here(), "/processed data/breast sample swapped ids in the first sequencing run_CCL.csv"))

path_raw <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "sequential_samples_hossein")
# cbc <- 
#   readxl::read_xlsx(paste0(path_raw, "/raw_data/breast_sarcama_report_Yifen_April2023.xlsx"),
#                     sheet = "CBC") %>% 
#   janitor::clean_names()
Demographic <- 
  readxl::read_xlsx(paste0(path_raw, "/raw_data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "Demographics", na = c("Missing", "Unknown")) %>% 
  janitor::clean_names()
# neutrophlil <- 
#   read_csv(paste0(path_save, "/processed data/Neutrophil lab data.csv"))

chart_dat <- 
  readxl::read_xlsx(paste0(#path_raw, 
    here::here(),
    "/chart_reviewed/Breast_chart reviews template_11252024.xlsx"),
    sheet = "Treatment_CBCs", na = c("NA", "UNK", "ND")) %>% 
  janitor::clean_names()

############################################################ II ### Merge data and limit to patients sequenced----
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
  mutate(sequenced_samples = "Yes") %>% 
  # add the swapped sample ids
  full_join(., fixed_ids %>% 
              mutate(`Sample Name (Secondary)` = received_ID), 
            by = c("MRN", "Sample Name (Secondary)", 
                   "Sample_Name_Fastq_file" = "Sample_Name")) %>% 
  select(MRN, Sample_Name_Fastq_file,
         sample_name_secondary = "Sample Name (Secondary)", 
         submitted_ID_for_added_samples = submitted_ID, 
         received_ID_for_added_samples = received_ID,
         status, treatment, sequenced_samples) %>% 
  # Fill the treatment for the extra samples that Nancy choose later
  group_by(MRN) %>% 
  fill(treatment, .direction = "updown") %>% 
  ungroup() %>% 
  distinct(MRN, Sample_Name_Fastq_file, .keep_all = TRUE)

sequenced_patients <- blood_patients %>% 
  # mutate(has_clinical = "Yes") %>% 
  full_join(., CH_status %>% 
              # mutate(has_CH = "Yes") %>% 
              mutate(sample_id = 
                       coalesce(submitted_ID_for_added_samples, sample_name_secondary)) %>% 
              mutate(sample_id = str_to_lower(sample_id)),
            by = c("mrn" = "MRN", "sample_id"
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
         # received_gcsf, gcsf_type, 
         neutropenia_at_anytime,
         ssf1_nm : progesterone_receptor_total_allred_score_desc,
         ER_results, PR_results, HER_results, 
         ER_PR_status, ER_PR_HER_status
  ) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  left_join(sequenced_patients, ., by = c("mrn", "date_of_birth" = "birth_dt")) %>% 
  mutate(race = coalesce(race, race_cancer_registry_1)) %>% 
  mutate(ethnicity = coalesce(ethnicity, ethnicity_cerner)) %>% 
  select(-age_at_diagnosis)

# Save intermediary data
# write_rds(sequenced_patients, paste0("sequenced_patients_", today(), ".rds"))

# rm
rm(blood_patients, later_blood_patients,
   submitted_samples, received_samples,
   CH_status, fixed_ids, fixed_ids_1,
   Demographic
   )


############################################################ III ### Update treatment dates and survival with chart review sept 2024----
# Cleaning chart data----
chart_dat <- chart_dat %>% 
  # Remove duplicated rows which appears between my output file and what Danny received
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
    # Set end date at 1 month from chart review for hormone end date in case of current use
    # Important to exclude CBC and other data
    hormone_end_date_current == "Currently use"     ~ as.Date("2024-09-23") + months(1) 
  )) %>% 
  select(-c(x0, x63))

# Sort samples vs new treatment dates
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
sequenced_patient_data <- sequenced_patients %>% 
  # mutate(specimen_collection_date = as.POSIXct(specimen_collection_date, origin = "1899-12-30")) %>% 
  mutate(sample_id = str_to_upper(sample_id)) %>% 
  select(-boost_dose_c_gy_1) %>% 
  full_join(.,
            chart_dat1, 
            by = c("mrn", 
                   "sample_name_secondary")) %>% 
  full_join(.,
            chart_dat %>% 
              select(mrn, type_of_radiation : boost_dose_c_gy_1) %>% 
              distinct(), 
            by = c("mrn")) %>% 
  # Add chart reviewed treatment clinical other than dates
  mutate(type_of_radiation_1 = type_of_radiation, .after = radiation_start_date_1) %>% 
  mutate(total_dose_radiation_1 = total_dose, .after = type_of_radiation_1) %>% 
  mutate(boost_dose_c_gy_radiation_1 = boost_dose_c_gy_1, .after = total_dose_radiation_1) %>% 
  select(mrn : Sample_Name_Fastq_file, party_id,
         specimen_collection_date, treatment_type,
         time_to_treatment, 
         sequence_sample_vs_treatment, treatment_received,
         everything(), 
         -c(type_of_radiation, total_dose, 
            boost_dose_c_gy_1, collection_date)) %>% 
  # create a discrepancies variables between patient's treatment type group
  # in case of new treatment were recorded
  mutate(discrepancies_old_new_sample_sequence = case_when(
    treatment_type == "chemo" &
      str_detect(sequence_sample_vs_treatment, "after_hormone|after_radiation")      ~ "need to be reclassified/remove",
    treatment_type == "chemorad" &
      str_detect(sequence_sample_vs_treatment, "after_hormone")                      ~ "need to be reclassified/remove",
    
  ), .after = treatment_type)


# END re-classified samples and update clinical dates with chart review ----

# Update survival and other clinical with chart review----
sequenced_patient_data <- chart_dat %>% 
  # Select 1 clinical row for each patient
  select(mrn, chemo1_end_date1_1,
         first_treatment_h_end_date,
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
              rename(date_of_diagnosis = date_of_diagnosis1) %>% 
              # Remove initially created variables before to add updated chart reviewed var
              select(-c(
                smoking_status, new_vital_status,
                contains("progesterone"), contains("estrogen"),
                type_of_first_recurrence_desc,
                ends_with("_results"), starts_with("her2"), 
                ER_PR_status, ER_PR_HER_status, 
                contains("oncotype_dx")
              )), 
            ., 
            by = "mrn") %>% 
  mutate(chart_reviewed_chemotherapy_end_date_1 = chemo1_end_date1_1,
         .after = chemotherapy_end_date1_1) %>% 
  mutate(chart_reviewed_end_date_of_first_corresponding_treatment = first_treatment_h_end_date,
         .after = chemo1_end_date1_1) %>% 
  select(-c(first_treatment_h_end_date, chemo1_end_date1_1,
            chemotherapy_end_date1_1))


# Save re-classified sequenced sequential sample and chart review clinical data
write_rds(sequenced_patient_data, 
          paste0(here::here(), 
                 "/processed data",
                 "/Identified breast data with re-classified sequenced sequential sample and chart review primary clinical_",
                 today(), ".rds"))


# END chart review data update
