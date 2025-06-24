library(tidyverse)
library(lubridate)

############################################################ I ### Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read_csv(paste0(
  here::here(),
  "/processed data",
  "/Identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia-NADIR_2025-05-16.csv"))

path_raw <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "sequential_samples_hossein")
pfs_chart_rev <- 
  readxl::read_xlsx(
    paste0(path_raw, "/chart_reviewed",
           "/updated PFS date chart review 06-24-25.xlsx")) %>% 
  janitor::clean_names()


############################################################ II ### Prep new data date----
pfs_chart_rev <- pfs_chart_rev %>% 
  `colnames<-`(c("deidentified_patient_id", "mrn",
                 "pfs_time_months_new", "date_of_progression_new", "pfs_date_new")) %>% 
  left_join(., 
            sequenced_patient_data %>% select(mrn, date_of_first_corresponding_treatment),
            by = "mrn") %>% 
  mutate(pfs_time_months_new = interval(start = date_of_first_corresponding_treatment, end = pfs_date_new)/
           duration(n = 1, unit = "months")) %>% 
  select(-date_of_first_corresponding_treatment) %>% 
  distinct()
  
  
############################################################ III ### Replace old date and length/interval variable----
sequenced_patient_data <- sequenced_patient_data %>% 
  left_join(., pfs_chart_rev,
            by = c("mrn", "deidentified_patient_id")) %>% 
  mutate(pfs_time_months = coalesce(pfs_time_months_new, pfs_time_months)) %>% 
  mutate(date_of_progression = coalesce(date_of_progression_new, date_of_progression)) %>% 
  mutate(pfs_date = coalesce(pfs_date_new, pfs_date)) %>% 
  select(-c(pfs_time_months_new, date_of_progression_new, pfs_date_new))

write_csv(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data_",
                                         today(), ".csv"))
write_csv(sequenced_patient_data, paste0(path_save, "/processed data",
                                         "/Identified breast data_",
                                         today(), ".csv"))
write_rds(sequenced_patient_data, paste0(path_save, "/processed data",
                                         "/Identified breast data_",
                                         today(), ".rds"))
write_rds(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data_",
                                         today(), ".rds"))


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
                 "/De-identified breast_", 
                 today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data", 
                 "/De-identified breast_", 
                 today(), ".csv"))

  
  
  
  
  
  
  
  
  