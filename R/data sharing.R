# Data sharing

library(tidyverse)

path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read_rds(paste0(# path_save,
  here::here(),
  "/processed data",
  "/Identified breast data_2025-06-30.rds"))

path_raw <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "CH_calls")

mutation_data <- readxl::read_xlsx(paste0(path_raw,
  "/Breast_final_calls_removeSynonymous_0827_08.12.2025.xlsx"), 
  sheet = "Sheet1")

survival_data <- mutation_data %>% 
  select(MRN, `included in Mona’s manuscript`) %>% 
  distinct() %>% 
  left_join(., sequenced_patient_data %>% 
              select(mrn, deidentified_patient_id,
                     date_of_first_corresponding_treatment,
                     date_of_last_followup_include_death_date = date_of_last_followup,
                     os_event, os_time_from_tx_start_months, 
                     pfs_date,
                     pfs_event, pfs_time_months) %>% 
              distinct(mrn, .keep_all = TRUE), 
            by = c("MRN" = "mrn"))

write_csv(survival_data, 
          paste0(path_save, 
                 "/processed data", 
                 "/Survival data matched with sample IDs from mutation data_", 
                 today(), ".csv"))


# End






