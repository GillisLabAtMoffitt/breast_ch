Global_data <- read_rds("/Users/colinccm/Documents/GitHub/Gillis/breast_ch/Global_data.rds")

blood_patients <- Global_data %>% 
  # filter to patients who have blood samples
  filter(!is.na(specimen_collection_date))

# blood_patients <- blood_patients %>% 
  # create age
  # mutate(age_at_diagnosis = interval(start = date_of_birth, end = date_of_diagnosis)/
  #          duration(n = 1, units = "years")) %>% 
  # create treatments cat
  # mutate(had_treatment = case_when(
  #   !is.na(chemotherapy_start_date_1) |
  #     !is.na(hormone_therapy_start_date_1) |
  #     !is.na(immunotherapy_start_date_1) |
  #     !is.na(radiation_start_date_1)                                ~ "Yes",
  #   TRUE                                                            ~ "No"
  # )) %>% 
  # mutate(blood_bf_chemo = case_when(
  #   specimen_collection_date <= chemotherapy_start_date_1           ~ "Yes",
  #   TRUE                                                            ~ "No"
  # )) %>% 
  # mutate(blood_bf_horm = case_when(
  #   specimen_collection_date <= hormone_therapy_start_date_1        ~ "Yes",
  #   TRUE                                                            ~ "No"
  # )) %>% 
  # mutate(blood_bf_immu = case_when(
  #   specimen_collection_date <= immunotherapy_start_date_1          ~ "Yes",
  #   TRUE                                                            ~ "No"
  # )) %>% 
  # mutate(blood_bf_rad = case_when(
  #   specimen_collection_date <= radiation_start_date_1              ~ "Yes",
  #   TRUE                                                            ~ "No"
  # )) %>% 
  # mutate(blood_bf_treatment = case_when(
  #   specimen_collection_date <= chemotherapy_start_date_1 &
  #     specimen_collection_date <= hormone_therapy_start_date_1 &
  #     specimen_collection_date <= immunotherapy_start_date_1 &
  #     specimen_collection_date <= radiation_start_date_1            ~ "Yes",
  #   TRUE                                                            ~ "No"
  # )) %>% 
  # mutate(blood_bf_treatment_30days = case_when(
  #   specimen_collection_date <= (chemotherapy_start_date_1 + days(30)) &
  #     specimen_collection_date <= (hormone_therapy_start_date_1 + days(30)) &
  #     specimen_collection_date <= (immunotherapy_start_date_1 + days(30)) &
  #     specimen_collection_date <= (radiation_start_date_1 + days(30))
  #   ~ "Yes",
  #   TRUE                                                            ~ "No"
  # )) %>% 
  
  
write_rds(blood_patients, "blood_patients.rds")


# End create variables
