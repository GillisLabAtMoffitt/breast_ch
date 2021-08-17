germline_data <- Global_data %>% 
  filter(!is.na(specimen_collection_date)) %>% 
  mutate(had_treatment = case_when(
    !is.na(chemotherapy_start_date_1) |
      !is.na(hormone_therapy_start_date_1) |
      !is.na(immunotherapy_start_date_1) |
      !is.na(radiation_start_date_1)                                ~ "Yes",
    TRUE                                                            ~ "No"
  )) %>% 
  mutate(blood_bf_chemo = case_when(
    specimen_collection_date <= chemotherapy_start_date_1           ~ "Yes",
    TRUE                                                            ~ "No"
  )) %>% 
  mutate(blood_bf_horm = case_when(
    specimen_collection_date <= hormone_therapy_start_date_1        ~ "Yes",
    TRUE                                                            ~ "No"
  )) %>% 
  mutate(blood_bf_immu = case_when(
    specimen_collection_date <= immunotherapy_start_date_1          ~ "Yes",
    TRUE                                                            ~ "No"
  )) %>% 
  mutate(blood_bf_rad = case_when(
    specimen_collection_date <= radiation_start_date_1              ~ "Yes",
    TRUE                                                            ~ "No"
  )) %>% 
  mutate(blood_bf_treatment = case_when(
    specimen_collection_date <= chemotherapy_start_date_1 &
      specimen_collection_date <= hormone_therapy_start_date_1 &
      specimen_collection_date <= immunotherapy_start_date_1 &
      specimen_collection_date <= radiation_start_date_1            ~ "Yes",
    TRUE                                                            ~ "No"
  )) #%>% 
  # mutate(blood_bf_treatment_30days = case_when(
  #   specimen_collection_date <= (chemotherapy_start_date_1 + 30) &
  #     specimen_collection_date <= (hormone_therapy_start_date_1 + 30) &
  #     specimen_collection_date <= (immunotherapy_start_date_1 + 30) &
  #     specimen_collection_date <= (radiation_start_date_1 + 30)       ~ "Yes",
  #   TRUE                                                            ~ "No"
  # ))




















write_rds(germline_data, "germline_data.rds")
