# Data sharing

library(tidyverse)

path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read_rds(paste0(# path_save,
  here::here(),
  "/processed data",
  "/Identified breast data_2025-08-15.rds"))


survival_data <- sequenced_patient_data %>% 
              select(mrn, deidentified_patient_id, 
                     patient_included_in_mona_manuscript, samples_excluded,
                     treatment_type,
                     date_of_first_corresponding_treatment,
                     date_of_last_followup_include_death_date = date_of_last_followup,
                     os_event, os_time_from_tx_start_months, 
                     pfs_date,
                     pfs_event, pfs_time_months) %>% 
              distinct(mrn, .keep_all = TRUE)

write_csv(survival_data, 
          paste0(path_save, 
                 "/processed data", 
                 "/Survival data matched with sample IDs from mutation data_", 
                 today(), ".csv"))


# End






