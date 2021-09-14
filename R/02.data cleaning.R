################################################################################# II ### Data cleaning
# Demo
Demographic <- Demographic %>% 
  mutate(across(where(is.character), ~str_to_lower(.)))


# Chemot
Chemot <- Chemot %>% 
  filter(chemotherapy_type == "CHEMO NOS" |
           chemotherapy_type == "MULTI-AGENT CHEMO" |
           chemotherapy_type == "NONE, NOT PLANNED" | # clean more
           chemotherapy_type == "RECOMMENDED,UNKN IF GIVEN" | # clean more
           chemotherapy_type == "SINGLE-AGENT CHEMO" |
           chemotherapy_type == "UNKNOWN; DC ONLY"  # clean more
  ) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # Make it easier to not have na for future filtering, rescue the ones with a drug name
  mutate(chemotherapy_completion_status_first_course = case_when(
    !is.na(chemotherapy_drug) |
      !is.na(chemotherapy_start_date)                        ~ coalesce(chemotherapy_completion_status_first_course, "chemo given")
  )) %>% 
  
  # filtering
  # Remove NA in both drug name and date
  filter(!(is.na(chemotherapy_drug) & is.na(chemotherapy_start_date) & is.na(chemotherapy_end_date))) %>%  # remove 341
  # Remove no chemo given in chemotherapy_drug
  filter(chemotherapy_drug != "no chemo given" | # remove 137
           is.na(chemotherapy_drug)) %>% 
  # Remove no chemotherapy when 18.. or 2300 dates only but keep if real date
  filter(!(str_detect(chemotherapy_start_date, "12:00:00 am") & # remove 271
             chemotherapy_completion_status_first_course == "no chemotherapy")) %>% 
  # To help after removing the unknown date
  # mutate(drugs_unk = case_when(
  #   str_detect(chemotherapy_start_date, "12:00:00 am") &
  #     is.na(chemotherapy_drug)                              ~ "unk drug, 1800 date",
  #   !is.na(chemotherapy_start_date) &
  #     is.na(chemotherapy_drug)                              ~ "unk drug",
  #   is.na(chemotherapy_start_date) &
  #     is.na(chemotherapy_drug)                              ~ "unk drug, NA date"
  # )) %>% 
  mutate(data_unk = case_when(
    str_detect(chemotherapy_start_date, "12:00:00 am") &
      !is.na(chemotherapy_drug)                             ~ "known drug, 1800 date",
    !is.na(chemotherapy_start_date) &
      !is.na(chemotherapy_drug)                             ~ "known drug",
    is.na(chemotherapy_start_date) &
      !is.na(chemotherapy_drug)                             ~ "unk drug, NA date",
    str_detect(chemotherapy_start_date, "12:00:00 am") &
      is.na(chemotherapy_drug)                              ~ "unk drug, 1800 date",
    !is.na(chemotherapy_start_date) &
      is.na(chemotherapy_drug)                              ~ "unk drug, known date",
  )) %>% 
  # Dates
  # remove the data 18.., 23..
  mutate(chemotherapy_start_date = case_when(
    str_detect(chemotherapy_start_date, "12:00:00 am")       ~ NA_character_,
    TRUE                                                     ~ chemotherapy_start_date
  ), 
  chemotherapy_start_date = as.Date(as.numeric(chemotherapy_start_date), 
                                    origin = "1899-12-30")
  ) %>% 
  mutate(chemotherapy_end_date = case_when(
    str_detect(chemotherapy_end_date, "12:00:00 am")         ~ NA_character_,
    TRUE                                                     ~ chemotherapy_end_date
  ), 
  chemotherapy_end_date = as.Date(as.numeric(chemotherapy_end_date), 
                                  origin = "1899-12-30")
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(chemotherapy_start_date = 
           coalesce(chemotherapy_start_date, 
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%  # 4,293 dates created
  # Fix the 2300 dates
  mutate(treatment_unk = case_when(
    str_detect(chemotherapy_start_date, "2300|2301")          ~ "Unknown when and if given 2300",
    TRUE                                                      ~ NA_character_
  )) %>% 
  mutate(chemotherapy_drug = coalesce(treatment_unk, chemotherapy_drug, data_unk)) %>% 
  mutate(chemotherapy_start_date = case_when(
    str_detect(chemotherapy_start_date, "2300|2301")          ~ as.Date("1700-01-01", origin = "1899-12-30"),
    TRUE                                                      ~ chemotherapy_start_date
  )) %>% 
  
  # Start pivot
  select(deidentified_patient_id, chemotherapy_drug, chemotherapy_start_date, chemotherapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, chemotherapy_start_date, chemotherapy_end_date) %>%
  summarise_at(vars(chemotherapy_drug), str_c, collapse = "; ") %>% 
  ungroup()

# pivot wider, use dcast bc better to keep date class
chemot <- dcast(setDT(Chemot), deidentified_patient_id ~ rowid(deidentified_patient_id),
                value.var = c(
                  "chemotherapy_drug",
                  "chemotherapy_start_date",
                  "chemotherapy_end_date"
                )
)

Chemot <- Chemot %>% 
  select(deidentified_patient_id, treatment_start_date = chemotherapy_start_date, 
         treatment_end_date = chemotherapy_end_date, treatment = chemotherapy_drug) %>% 
  mutate(treatment_type = "chemo")



# Hormonet
Hormonet <- Hormonet %>% 
  filter(hormone_therapy_type == "HORMONE ADMINISTERED" |
           hormone_therapy_type == "NONE, NOT PLANNED" |
           hormone_therapy_type == "RECOMMENDED,UNKN IF GIVEN" | # clean more
           hormone_therapy_type == "UNKNOWN; DC ONLY"  # clean more
  ) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  filter(hormone_therapy_drug != "no hormone given" | # remove 99
           is.na(hormone_therapy_drug)) %>% 
  
  mutate(data_unk = case_when(
    str_detect(hormone_therapy_start_date, "12:00:00 am") &
      !is.na(hormone_therapy_drug)                             ~ "known drug, 1800 date",
    !is.na(hormone_therapy_start_date) &
      !is.na(hormone_therapy_drug)                             ~ "known drug",
    is.na(hormone_therapy_start_date) &
      !is.na(hormone_therapy_drug)                             ~ "unk drug, NA date",
    str_detect(hormone_therapy_start_date, "12:00:00 am") &
      is.na(hormone_therapy_drug)                              ~ "unk drug, 1800 date",
    !is.na(hormone_therapy_start_date) &
      is.na(hormone_therapy_drug)                              ~ "unk drug, known date",
  )) %>% 
  filter(!is.na(data_unk)) %>% 
  # Dates
  # remove the data 18.., 23..
  mutate(hormone_therapy_start_date = case_when(
    str_detect(hormone_therapy_start_date, "12:00:00 am")       ~ NA_character_,
    TRUE                                                     ~ hormone_therapy_start_date
  ), 
  hormone_therapy_start_date = as.Date(as.numeric(hormone_therapy_start_date), 
                                       origin = "1899-12-30")
  ) %>% 
  mutate(hormone_therapy_end_date = case_when(
    str_detect(hormone_therapy_end_date, "12:00:00 am")         ~ NA_character_,
    TRUE                                                     ~ hormone_therapy_end_date
  ), 
  hormone_therapy_end_date = as.Date(as.numeric(hormone_therapy_end_date), 
                                     origin = "1899-12-30")
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(hormone_therapy_start_date = 
           coalesce(hormone_therapy_start_date, 
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # Fix the 2300 dates
  mutate(treatment_unk = case_when(
    str_detect(hormone_therapy_start_date, "2300|2301")          ~ "Unknown when and if given 2300",
    TRUE                                                      ~ NA_character_
  )) %>% 
  mutate(hormone_therapy_drug = coalesce(treatment_unk, hormone_therapy_drug, data_unk)) %>% 
  mutate(hormone_therapy_start_date = case_when(
    str_detect(hormone_therapy_start_date, "2300|2301")          ~ as.Date("1700-01-01", origin = "1899-12-30"),
    TRUE                                                      ~ hormone_therapy_start_date
  )) %>% 
  
  # Start pivot
  select(deidentified_patient_id, hormone_therapy_drug, hormone_therapy_start_date, hormone_therapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, hormone_therapy_start_date, hormone_therapy_end_date) %>%
  summarise_at(vars(hormone_therapy_drug), str_c, collapse = "; ") %>% 
  ungroup()

hormonet <- dcast(setDT(Hormonet), deidentified_patient_id ~ rowid(deidentified_patient_id),
                  value.var = c(
                    "hormone_therapy_drug",
                    "hormone_therapy_start_date",
                    "hormone_therapy_end_date"
                  )
)

Hormonet <- Hormonet %>% 
  select(deidentified_patient_id, treatment_start_date = hormone_therapy_start_date, 
         treatment_end_date = hormone_therapy_end_date, treatment = hormone_therapy_drug) %>% 
  mutate(treatment_type = "hormone")



# Immnunot
Immnunot <- Immnunot %>% 
  filter(immunotherapy_type == "IMMUNO ADMINISTERED" |
           immunotherapy_type == "NONE, NOT PLANNED" |
           immunotherapy_type == "RECOMMENDED,UNKN IF GIVEN" | # clean more
           immunotherapy_type == "UNKNOWN; DC ONLY"  # clean more
  ) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  mutate(immunotherapy_drug = case_when(
    is.na(immunotherapy_start_date) &
      !is.na(immunotherapy_drug)                             ~ "unk drug, NA date",
    str_detect(immunotherapy_start_date, "12:00:00 am") &
      is.na(immunotherapy_drug)                              ~ "unk drug, 1800 date",
    !is.na(immunotherapy_start_date) &
      is.na(immunotherapy_drug)                              ~ "unk drug, known date",
  )) %>% 
  filter(!is.na(immunotherapy_drug)) %>% 
  # Dates
  # remove the data 18.., 23..
  mutate(immunotherapy_start_date = case_when(
    str_detect(immunotherapy_start_date, "12:00:00 am")       ~ NA_character_,
    TRUE                                                     ~ immunotherapy_start_date
  ), 
  immunotherapy_start_date = as.Date(as.numeric(immunotherapy_start_date), 
                                     origin = "1899-12-30")
  ) %>% 
  mutate(immunotherapy_end_date = case_when(
    str_detect(immunotherapy_end_date, "12:00:00 am")         ~ NA_character_,
    TRUE                                                     ~ immunotherapy_end_date
  ), 
  immunotherapy_end_date = as.Date(as.numeric(immunotherapy_end_date),
                                   origin = "1899-12-30")
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(immunotherapy_start_date = 
           coalesce(immunotherapy_start_date, 
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # No 2300 date
  
  # Start pivot
  select(deidentified_patient_id, immunotherapy_drug, immunotherapy_start_date, immunotherapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, immunotherapy_start_date, immunotherapy_end_date) %>%
  summarise_at(vars(immunotherapy_drug), str_c, collapse = "; ") %>% 
  ungroup()

immnunot <- dcast(setDT(Immnunot), deidentified_patient_id ~ rowid(deidentified_patient_id),
                  value.var = c(
                    "immunotherapy_drug",
                    "immunotherapy_start_date",
                    "immunotherapy_end_date"
                  )
)

Immnunot <- Immnunot %>% 
  select(deidentified_patient_id, treatment_start_date = immunotherapy_start_date, 
         treatment_end_date = immunotherapy_end_date, treatment = immunotherapy_drug) %>% 
  mutate(treatment_type = "immunoT")



# Radiot
Radiot <- Radiot %>% 
  filter(reason_for_no_radiation == "RAD THERAPY PERFORMED"| 
           reason_for_no_radiation == "RECOMMENDED,UNKN IF GIVEN" | 
           reason_for_no_radiation == "UNKNOWN; DC ONLY" |
           is.na(reason_for_no_radiation)) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  
  mutate(data_unk = case_when(
    str_detect(radiation_start_date, "12:00:00 am") &
      !is.na(boost_dose_c_gy)                             ~ "known dose, 1800 date",
    !is.na(radiation_start_date) &
      !is.na(boost_dose_c_gy)                             ~ "known dose",
    is.na(radiation_start_date) &
      !is.na(boost_dose_c_gy)                             ~ "unk dose, NA date",
    str_detect(radiation_start_date, "12:00:00 am") &
      is.na(boost_dose_c_gy)                              ~ "unk dose, 1800 date",
    !is.na(radiation_start_date) &
      is.na(boost_dose_c_gy)                              ~ "unk dose, known date",
  )) %>% 
  filter(!is.na(data_unk)) %>% 
  # Dates
  # remove the data 18.., 23..
  mutate(radiation_start_date = case_when(
    str_detect(radiation_start_date, "12:00:00 am")       ~ NA_character_,
    TRUE                                                     ~ radiation_start_date
  ), 
  radiation_start_date = as.Date(as.numeric(radiation_start_date), 
                                 origin = "1899-12-30")
  ) %>% 
  mutate(radiation_end_date = case_when(
    str_detect(radiation_end_date, "12:00:00 am")         ~ NA_character_,
    TRUE                                                     ~ radiation_end_date
  ), 
  radiation_end_date = as.Date(as.numeric(radiation_end_date), ################################################ Fix bug create a 2300 date, lubridate::origin
                               origin = "1899-12-30")
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(radiation_start_date = 
           coalesce(radiation_start_date, 
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # Fix the 2300 dates
  mutate(treatment_unk = case_when(
    str_detect(radiation_start_date, "2300|2301")             ~ "Unknown when and if given 2300",
    TRUE                                                      ~ NA_character_
  )) %>%
  mutate(boost_dose_c_gy = coalesce(treatment_unk, as.character(boost_dose_c_gy), data_unk)) %>% 
  mutate(radiation_start_date = case_when(
    str_detect(radiation_start_date, "2300|2301")             ~ as.Date("1700-01-01", origin = "1899-12-30"),
    TRUE                                                      ~ radiation_start_date
  )) %>%
  
  # Start pivot
  select(deidentified_patient_id, boost_dose_c_gy, radiation_start_date, radiation_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, radiation_start_date, radiation_end_date) %>%
  summarise_at(vars(boost_dose_c_gy), str_c, collapse = "; ") %>% 
  ungroup()


radiot <- dcast(setDT(Radiot), deidentified_patient_id ~ rowid(deidentified_patient_id),
                value.var = c(
                  "boost_dose_c_gy",
                  "radiation_start_date",
                  "radiation_end_date"
                )
)


Radiot <- Radiot %>% 
  select(deidentified_patient_id, treatment_start_date = radiation_start_date, 
         treatment_end_date = radiation_end_date, treatment = boost_dose_c_gy) %>% 
  mutate(treatment_type = "radioT")




# Combine
treatment <- bind_rows(Chemot, Hormonet, Immnunot, Radiot) %>% 
  arrange(deidentified_patient_id, treatment_start_date) %>% 
  group_by(deidentified_patient_id, treatment_type) %>% 
  mutate(treatment_line = row_number(deidentified_patient_id)) %>% 
  unite(treatment_line, c(treatment_type, treatment_line), sep = "_", remove = FALSE)

Treatment <- full_join(chemot, hormonet, by = "deidentified_patient_id") %>% 
  full_join(., immnunot, by = "deidentified_patient_id") %>% 
  full_join(., radiot, by = "deidentified_patient_id") %>% 
  mutate(had_chemo = case_when(
    !is.na(chemotherapy_start_date_1) ~ "Yes"
  )) %>% 
  mutate(had_hormone = case_when(
    !is.na(hormone_therapy_start_date_1) ~ "Yes"
  )) %>% 
  mutate(had_immuno = case_when(
    !is.na(immunotherapy_start_date_1) ~ "Yes"
  )) %>% 
  mutate(had_rad = case_when(
    !is.na(radiation_start_date_1) ~ "Yes"
  )) %>% 
  mutate(had_chemo_rad = case_when(
    !is.na(chemotherapy_start_date_1) &
      !is.na(radiation_start_date_1) ~ "Yes"
  )) %>% 
  mutate(had_treatment = case_when(
    !is.na(chemotherapy_start_date_1) &
      !is.na(hormone_therapy_start_date_1) &
      !is.na(immunotherapy_start_date_1) &
      !is.na(radiation_start_date_1)             ~ "Yes"
  ))

write_rds(Treatment, "Treatment.rds")



# DNA, clean and add same sample/same date on the same row
breast_dna <- breast_DNA %>% 
  filter(derived_tissue_type == "Blood",sample_type != "WBC/RBC") %>% 
  mutate(deidentified_patient_id = str_to_lower(deidentified_patient_id)) %>% 
  select(deidentified_patient_id, sample_family_id_sf, sample_id,
         specimen_collection_date) %>%
  mutate(specimen_collection_date = as.Date(specimen_collection_date, format = "%Y-%M-%d")) %>% 
  arrange(deidentified_patient_id, specimen_collection_date) %>% 
  group_by(deidentified_patient_id, sample_family_id_sf, specimen_collection_date) %>% 
  summarise_at(vars(sample_id), str_c, collapse = "; ") %>%
  # separate(col = sample_id, paste("sample_id", 1:3, sep="_"), sep = "; ", extra = "drop", fill = "right")
  ungroup() %>% 
  left_join(breast_DNA)

# write_rds(breast_dna, "breast_dna.rds")

breast_dna1 <- breast_dna %>% left_join(., Treatment, by = "deidentified_patient_id") %>% 
  mutate(blood_bf_chemo = case_when(
    specimen_collection_date <= chemotherapy_start_date_1                ~ "Yes",
    specimen_collection_date > chemotherapy_start_date_1                 ~ "No",
    is.na(chemotherapy_start_date_1)                                     ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>% 
  mutate(blood_bf_hormone = case_when(
    specimen_collection_date <= hormone_therapy_start_date_1                ~ "Yes",
    specimen_collection_date > hormone_therapy_start_date_1                 ~ "No",
    is.na(hormone_therapy_start_date_1)                                     ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>% 
  mutate(blood_bf_immuno = case_when(
    specimen_collection_date <= immunotherapy_start_date_1                ~ "Yes",
    specimen_collection_date > immunotherapy_start_date_1                 ~ "No",
    is.na(immunotherapy_start_date_1)                                     ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>% 
  mutate(blood_bf_rad = case_when(
    specimen_collection_date <= radiation_start_date_1                ~ "Yes",
    specimen_collection_date > radiation_start_date_1                 ~ "No",
    is.na(radiation_start_date_1)                                     ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>% 
  mutate(blood_bf_chemo_rad = case_when(
    specimen_collection_date <= chemotherapy_start_date_1 &
      specimen_collection_date <= radiation_start_date_1                ~ "Yes",
    specimen_collection_date > chemotherapy_start_date_1 |
      specimen_collection_date > radiation_start_date_1                 ~ "No",
    is.na(chemotherapy_start_date_1) |
      is.na(radiation_start_date_1)                                 ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>% 
  mutate(blood_bf_treatment = case_when(
    specimen_collection_date <= chemotherapy_start_date_1 &
      specimen_collection_date <= hormone_therapy_start_date_1 &
      specimen_collection_date <= immunotherapy_start_date_1 &
      specimen_collection_date <= radiation_start_date_1                ~ "Yes",
    if_any(contains("start_date_1"), ~ . < specimen_collection_date)    ~ "No",
    # if_all(contains("start_date_1"), ~ . < specimen_collection_date)    ~ "Nope",
    is.na(chemotherapy_start_date_1) |
      is.na(hormone_therapy_start_date_1) |
      is.na(immunotherapy_start_date_1) |
      is.na(radiation_start_date_1)                                     ~ "not administred",
    TRUE                                                                ~ NA_character_
  )) %>% 
  mutate(blood_bf_30_days_chemo = case_when(
    specimen_collection_date >= (chemotherapy_start_date_1 - days(30)) &
      specimen_collection_date <= (chemotherapy_start_date_1 + days(30))              ~ "Yes",
    specimen_collection_date > (chemotherapy_start_date_1 + days(30))                 ~ "No",
    is.na(chemotherapy_start_date_1)                                     ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>%
  mutate(blood_bf_30_days_hormone = case_when(
    specimen_collection_date <= (hormone_therapy_start_date_1 + days(30))                ~ "Yes",
    specimen_collection_date > (hormone_therapy_start_date_1 + days(30))                 ~ "No",
    is.na(hormone_therapy_start_date_1)                                     ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>%
  mutate(blood_bf_30_days_immuno = case_when(
    specimen_collection_date <= (immunotherapy_start_date_1 + days(30))                ~ "Yes",
    specimen_collection_date > (immunotherapy_start_date_1 + days(30))                 ~ "No",
    is.na(immunotherapy_start_date_1)                                     ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>%
  mutate(blood_bf_30_days_rad = case_when(
    specimen_collection_date <= (radiation_start_date_1 + days(30))                ~ "Yes",
    specimen_collection_date > (radiation_start_date_1 + days(30))                 ~ "No",
    is.na(radiation_start_date_1)                                     ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>%
  mutate(blood_bf_30_days_chemo_rad = case_when(
    specimen_collection_date <= (chemotherapy_start_date_1 + days(30)) &
      specimen_collection_date <= (radiation_start_date_1 + days(30))                ~ "Yes",
    specimen_collection_date > (chemotherapy_start_date_1 + days(30)) &
      specimen_collection_date > (radiation_start_date_1 + days(30))                 ~ "No",
    is.na(chemotherapy_start_date_1) |
      is.na(radiation_start_date_1)                                 ~ "not administred",
    TRUE                                                            ~ NA_character_
  )) %>%
  mutate(blood_bf_30_days_treatment = case_when(
    specimen_collection_date <= (chemotherapy_start_date_1 + days(30)) &
      specimen_collection_date <= (hormone_therapy_start_date_1 + days(30)) &
      specimen_collection_date <= (immunotherapy_start_date_1 + days(30)) &
      specimen_collection_date <= (radiation_start_date_1 + days(30))                ~ "Yes",
    if_any(contains("start_date_1"), ~ (.  + days(30)) < specimen_collection_date)    ~ "No",
    # if_all(contains("start_date_1"), ~ . < specimen_collection_date)    ~ "Nope",
    is.na(chemotherapy_start_date_1) |
      is.na(hormone_therapy_start_date_1) |
      is.na(immunotherapy_start_date_1) |
      is.na(radiation_start_date_1)                                     ~ "not administred",
    TRUE                                                                ~ NA_character_
  )) %>%
  mutate(across(contains("blood_bf_"), ~ factor(., levels = c("Yes", "No", "not administred")))) %>% 
  
  # ungroup() %>% 
  # mutate(treatment_bf_blood = case_when(
  #   specimen_collection_date <= treatment_start_date                ~ "Blood first",
  #   specimen_collection_date <= treatment_start_date                ~ treatment_line,
  #   TRUE                                                            ~ NA_character_
  # )) %>% 
  # mutate(treatment_bf_30days_blood = case_when(
  #   specimen_collection_date <= (treatment_start_date + days(30))   ~ "Blood first",
  #   specimen_collection_date <= (treatment_start_date + days(30))   ~ treatment_line,
  #   TRUE                                                            ~ NA_character_
# )) %>% 
# 
# mutate(treatment_after_blood = case_when(
#   specimen_collection_date > treatment_start_date                 ~ treatment_line,
#   TRUE                                                            ~ NA_character_
# )) %>% 
# mutate(treatment_after_30days_blood = case_when(
#   specimen_collection_date > (treatment_start_date + days(30))   ~ treatment_line,
#   TRUE                                                            ~ NA_character_
# )) %>% 
# distinct(deidentified_patient_id, sample_family_id_sf, sample_id,
#          specimen_collection_date, )
arrange(deidentified_patient_id, specimen_collection_date, blood_bf_treatment) %>% 
  group_by(deidentified_patient_id) %>% 
  # mutate(sample_lag = specimen_collection_date - lag(specimen_collection_date), 
  #        sample_lag = str_remove(sample_lag, " days")
  #        ) %>% 
  # mutate(has_a_good_sample = case_when(
  #   if_any(contains("blood_bf_"), ~ . == "Yes")    ~ "Yes"
  # )) %>% 
  # fill(has_a_good_sample, .direction = "down") %>% 
  # mutate(has_a_good_seq_sample = case_when(
  #   blood_bf_30_days_chemo_rad == "No" &
  #     has_a_good_sample == "Yes" &
  #     sample_lag > 30                             ~ "Yes"
# )) %>% 

select(deidentified_patient_id, specimen_collection_date, #sample_lag, has_a_good_sample, has_a_good_seq_sample,
       "blood_bf_chemo", "blood_bf_hormone", 
       "blood_bf_immuno", "blood_bf_rad",
       "blood_bf_chemo_rad",
       "blood_bf_treatment", everything())

# breast_dna1 %>% distinct(deidentified_patient_id, .keep_all = TRUE) %>% 
#   select(blood_bf_chemo, blood_bf_hormone, 
#          blood_bf_immuno, blood_bf_rad,
#          blood_bf_treatment, blood_bf_30_days_treatment) %>% 
#   tbl_summary(sort = list(everything() ~ "frequency"))

breast_dna1 %>% filter(chemotherapy_start_date_1 == "1700-01-01") %>% nrow()
breast_dna1 %>% filter(hormone_therapy_start_date_1 == "1700-01-01") %>% nrow()
breast_dna1 %>% filter(immunotherapy_start_date_1 == "1700-01-01") %>% nrow()
breast_dna1 %>% filter(radiation_start_date_1 == "1700-01-01") %>% nrow()
breast_dna1 %>% filter(chemotherapy_start_date_1 == "1700-01-01" & 
                         radiation_start_date_1 == "1700-01-01") %>% nrow()
breast_dna1 %>% filter(chemotherapy_start_date_1 == "1700-01-01" & 
                         hormone_therapy_start_date_1 == "1700-01-01" &
                         immunotherapy_start_date_1 == "1700-01-01" &
                         radiation_start_date_1 == "1700-01-01") %>% nrow()



# For samples before treatment, take only the closest

sample_before_chemo <- breast_dna1 %>% 
  filter(blood_bf_chemo == "Yes") %>% 
  # pivot_longer(cols = c(chemotherapy_start_date_1, hormone_therapy_start_date_1, 
  #                       immunotherapy_start_date_1, radiation_start_date_1),
  #              names_to = "treatmetn_type", values_to = "date_treatment") %>% 
  mutate(interval_sample_chemo = 
           abs(interval(start = specimen_collection_date, end = chemotherapy_start_date_1) /
                 duration(n = 1, units = "days"))) %>% 
  arrange(deidentified_patient_id, interval_sample_chemo) %>% 
  distinct(deidentified_patient_id, .keep_all = TRUE)

library(ggforce)
p <- qplot(x =interval_sample_chemo, data=subset(sample_before_chemo), fill=..count.., 
           geom="histogram", 
           binwidth = 100,
) 
p + scale_fill_viridis_c(
  alpha = 1,
  begin = 0,
  end = 1,
  direction = 1,
  option = "D",
  values = NULL,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
) +
  theme_minimal(base_size = 14) +
  labs(x="Time from Blood Collection to Chemotherapy (in days)", 
       y="Number of Patient",
       caption = "Each bar represents 100 days") +
  facet_zoom(ylim = c(0, 10), zoom.size = 1)

sample_before_hormone <- breast_dna1 %>% 
  filter(blood_bf_hormone == "Yes") %>% 
  mutate(interval_sample_hormone = 
           abs(interval(start = specimen_collection_date, end = hormone_therapy_start_date_1) /
                 duration(n = 1, units = "days"))) %>% 
  arrange(deidentified_patient_id, interval_sample_hormone) %>% 
  distinct(deidentified_patient_id, .keep_all = TRUE)

p <- qplot(x =interval_sample_hormone, data=subset(sample_before_hormone), fill=..count.., 
           geom="histogram", 
           binwidth = 100,
) 
p + scale_fill_viridis_c(
  alpha = 1,
  begin = 0,
  end = 1,
  direction = 1,
  option = "A",
  values = NULL,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
) +
  theme_minimal(base_size = 14) +
  labs(x="Time from Blood Collection to Hormonetherapy (in days)", 
       y="Number of Patient",
       caption = "Each bar represents 100 days") +
  facet_zoom(ylim = c(0, 10), zoom.size = 1)

sample_before_immuno <- breast_dna1 %>% 
  filter(blood_bf_immuno == "Yes") %>% 
  mutate(interval_sample_immuno = 
           abs(interval(start = specimen_collection_date, end = immunotherapy_start_date_1) /
                 duration(n = 1, units = "days"))) %>% 
  arrange(deidentified_patient_id, interval_sample_immuno) %>% 
  distinct(deidentified_patient_id, .keep_all = TRUE)

p <- qplot(x =interval_sample_immuno, data=subset(sample_before_immuno), fill=..count.., 
           geom="histogram", 
           binwidth = 100,
) 
p + scale_fill_viridis_c(
  alpha = 1,
  begin = 0,
  end = 1,
  direction = 1,
  option = "H",
  values = NULL,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
) +
  theme_minimal(base_size = 14) +
  labs(x="Time from Blood Collection to Immunotherapy (in days)", 
       y="Number of Patient",
       caption = "Each bar represents 100 days") +
  facet_zoom(ylim = c(0, 10), zoom.size = 1)

sample_before_rad <- breast_dna1 %>% 
  filter(blood_bf_rad == "Yes") %>% 
  mutate(interval_sample_rad = 
           abs(interval(start = specimen_collection_date, end = radiation_start_date_1) /
                 duration(n = 1, units = "days"))) %>% 
  arrange(deidentified_patient_id, interval_sample_rad) %>% 
  distinct(deidentified_patient_id, .keep_all = TRUE)

p <- qplot(x =interval_sample_rad, data=subset(sample_before_rad), fill=..count.., 
           geom="histogram", 
           binwidth = 100,
) 
p + scale_fill_viridis_c(
  alpha = 1,
  begin = 0,
  end = 1,
  direction = 1,
  option = "A",
  values = NULL,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
) +
  theme_minimal(base_size = 14) +
  labs(x="Time from Blood Collection to Radiotherapy (in days)", 
       y="Number of Patient",
       caption = "Each bar represents 100 days") +
  facet_zoom(ylim = c(0, 5), zoom.size = 1)

sample_before <- bind_rows(sample_before_chemo, sample_before_hormone, sample_before_immuno, sample_before_rad)



breast_dna2 <- breast_dna1 %>% 
  mutate(had_good_sample_chemo = case_when(
    blood_bf_chemo == "Yes"              ~ "Yes"
  )) %>% 
  group_by(deidentified_patient_id) %>% 
  fill(had_good_sample_chemo, .direction = "updown") %>% 
  mutate(seq_sample_chemo = case_when(
    had_good_sample_chemo == "Yes" &
      blood_bf_chemo == "No" ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  mutate(had_good_sample_hormone = case_when(
    blood_bf_hormone == "Yes"              ~ "Yes"
  )) %>% 
  group_by(deidentified_patient_id) %>% 
  fill(had_good_sample_hormone, .direction = "updown") %>% 
  mutate(seq_sample_hormone = case_when(
    had_good_sample_hormone == "Yes" &
      blood_bf_hormone == "No" ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  mutate(had_good_sample_immuno = case_when(
    blood_bf_immuno == "Yes"              ~ "Yes"
  )) %>% 
  group_by(deidentified_patient_id) %>% 
  fill(had_good_sample_immuno, .direction = "updown") %>% 
  mutate(seq_sample_immuno = case_when(
    had_good_sample_immuno == "Yes" &
      blood_bf_immuno == "No" ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  mutate(had_good_sample_rad = case_when(
    blood_bf_rad == "Yes"              ~ "Yes"
  )) %>% 
  group_by(deidentified_patient_id) %>% 
  fill(had_good_sample_rad, .direction = "updown") %>% 
  mutate(seq_sample_rad = case_when(
    had_good_sample_rad == "Yes" &
      blood_bf_rad == "No" ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  mutate(had_good_sample_chemo_rad = case_when(
    blood_bf_chemo_rad == "Yes"              ~ "Yes"
  )) %>% 
  group_by(deidentified_patient_id) %>% 
  fill(had_good_sample_chemo_rad, .direction = "updown") %>% 
  mutate(seq_sample_chemo_rad = case_when(
    had_good_sample_chemo_rad == "Yes" &
      blood_bf_chemo_rad == "No" ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  mutate(had_good_sample_treatment = case_when(
    blood_bf_treatment == "Yes"              ~ "Yes"
  )) %>% 
  group_by(deidentified_patient_id) %>% 
  fill(had_good_sample_treatment, .direction = "updown") %>% 
  mutate(seq_sample_treatment = case_when(
    had_good_sample_treatment == "Yes" &
      blood_bf_treatment == "No" ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  ungroup()




# left_join(., Treatment %>% 
#             select(deidentified_patient_id, had_chemo, had_hormone, had_immuno, had_rad, had_chemo_rad, had_treatment),
#           by = "deidentified_patient_id")

sample_after_chemo <- breast_dna2 %>%
  filter(seq_sample_chemo == "Yes") %>%
  mutate(interval_sample_chemo =
           abs(interval(start = specimen_collection_date, end = chemotherapy_start_date_1) /
                 duration(n = 1, units = "days"))) %>%
  arrange(deidentified_patient_id, interval_sample_chemo) %>%
  distinct(deidentified_patient_id, specimen_collection_date, .keep_all = TRUE) %>% 
  group_by(deidentified_patient_id) %>% 
  mutate(sequential_sample_count = factor(row_number(deidentified_patient_id))) %>% 
  ungroup() %>% 
  select(deidentified_patient_id, sequential_sample_count, interval_sample_chemo, everything())

sample_after_chemo %>% 
  ggplot(aes(x =interval_sample_chemo, fill = sequential_sample_count))+
  geom_histogram(binwidth = 100, alpha = 0.9, position = "stack") +
  scale_fill_viridis(discrete=T)+
  theme_minimal(base_size = 14) +
  labs(x="Time from Blood Collection to Chemotherapy (in days)", 
       y="Number of Patient",
       caption = "Each bar represents 100 days")

sample_after_hormone <- breast_dna2 %>%
  filter(seq_sample_hormone == "Yes") %>%
  mutate(interval_sample_hormone =
           abs(interval(start = specimen_collection_date, end = hormone_therapy_start_date_1) /
                 duration(n = 1, units = "days"))) %>%
  arrange(deidentified_patient_id, interval_sample_hormone) %>%
  distinct(deidentified_patient_id, specimen_collection_date, .keep_all = TRUE) %>% 
  group_by(deidentified_patient_id) %>% 
  mutate(sequential_sample_count = factor(row_number(deidentified_patient_id))) %>% 
  ungroup() %>% 
  select(deidentified_patient_id, sequential_sample_count, interval_sample_hormone, everything())

sample_after_hormone %>% 
  ggplot(aes(x =interval_sample_hormone, fill = sequential_sample_count))+
  geom_histogram(binwidth = 100, alpha = 0.9, position = "stack") +
  scale_fill_viridis(discrete=T)+
  theme_minimal(base_size = 14) +
  labs(x="Time from Blood Collection to Hormonetherapy (in days)", 
       y="Number of Patient",
       caption = "Each bar represents 100 days")

sample_after_immuno <- breast_dna2 %>%
  filter(seq_sample_immuno == "Yes") %>%
  mutate(interval_sample_immuno =
           abs(interval(start = specimen_collection_date, end = immunotherapy_start_date_1) /
                 duration(n = 1, units = "days"))) %>%
  arrange(deidentified_patient_id, interval_sample_immuno) %>%
  distinct(deidentified_patient_id, specimen_collection_date, .keep_all = TRUE) %>% 
  group_by(deidentified_patient_id) %>% 
  mutate(sequential_sample_count = factor(row_number(deidentified_patient_id))) %>% 
  ungroup() %>% 
  select(deidentified_patient_id, sequential_sample_count, interval_sample_immuno, everything())

sample_after_immuno %>% 
  ggplot(aes(x =interval_sample_immuno, fill = sequential_sample_count))+
  geom_histogram(binwidth = 100, alpha = 0.9, position = "stack") +
  scale_fill_viridis(discrete=T)+
  theme_minimal(base_size = 14) +
  labs(x="Time from Blood Collection to Immunotherapy (in days)", 
       y="Number of Patient",
       caption = "Each bar represents 100 days")

sample_after_rad <- breast_dna2 %>%
  filter(seq_sample_rad == "Yes") %>%
  mutate(interval_sample_rad =
           abs(interval(start = specimen_collection_date, end = radiation_start_date_1) /
                 duration(n = 1, units = "days"))) %>%
  arrange(deidentified_patient_id, interval_sample_rad) %>%
  distinct(deidentified_patient_id, specimen_collection_date, .keep_all = TRUE) %>% 
  group_by(deidentified_patient_id) %>% 
  mutate(sequential_sample_count = factor(row_number(deidentified_patient_id))) %>% 
  ungroup() %>% 
  select(deidentified_patient_id, sequential_sample_count, interval_sample_rad, everything())

sample_after_rad %>% 
  ggplot(aes(x =interval_sample_rad, fill = sequential_sample_count))+
  geom_histogram(binwidth = 100, alpha = 0.9, position = "stack") +
  scale_fill_viridis(discrete=T)+
  theme_minimal(base_size = 14) +
  labs(x="Time from Blood Collection to Radiotherapy (in days)", 
       y="Number of Patient",
       caption = "Each bar represents 100 days")

sample_after <- bind_rows(sample_after_chemo, sample_after_hormone, sample_after_immuno, sample_after_rad)













# breast_dna2 <-
#   dcast(setDT(breast_dna1), deidentified_patient_id ~ rowid(deidentified_patient_id),
#         value.var = c("sample_family_id_sf", 
#                       "sample_id",
#                       "specimen_collection_date", 
#                       "blood_bf_chemo", "blood_bf_hormone", 
#                       "blood_bf_immuno", "blood_bf_rad",
#                       "blood_bf_chemo_rad",
#                       "blood_bf_treatment", 
#                       "blood_bf_30_days_chemo", "blood_bf_30_days_hormone", 
#                       "blood_bf_30_days_immuno", "blood_bf_30_days_rad",
#                       "blood_bf_30_days_chemo_rad",
#                       "blood_bf_30_days_treatment"),
#         sep = "_sample") %>% 
#   
#   left_join(., Treatment %>% 
#               select(deidentified_patient_id, had_chemo, had_hormone, had_immuno, had_rad, had_chemo_rad, had_treatment),
#             by = "deidentified_patient_id") %>% 
#   
#   mutate(seq_sample_chemo = case_when(
#     blood_bf_chemo_sample1 == "Yes" &
#       if_any(contains("blood_bf_chemo_sample"), ~ . == "No") ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(seq_sample_hormone = case_when(
#     blood_bf_hormone_sample1 == "Yes" &
#       if_any(contains("blood_bf_hormone_sample"), ~ . == "No") ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(seq_sample_immuno = case_when(
#     blood_bf_immuno_sample1 == "Yes" &
#       if_any(contains("blood_bf_immuno_sample"), ~ . == "No") ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(seq_sample_rad = case_when(
#     blood_bf_rad_sample1 == "Yes" &
#       if_any(contains("blood_bf_rad_sample"), ~ . == "No") ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(seq_sample_chemo_rad = case_when(
#     blood_bf_chemo_rad_sample1 == "Yes" &
#       if_any(contains("blood_bf_chemo_rad_sample"), ~ . == "No") ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(seq_sample_treatment = case_when(
#     blood_bf_treatment_sample1 == "Yes" &
#       if_any(contains("blood_bf_treatment_sample"), ~ . == "No") ~ "Yes",
#     TRUE ~ "No"
#   )) 
# 
# 
# 
# 
# # Table if take samples strictly before treatment and sequential during 30 days after treatment
# tbl1 <- breast_dna2 %>% #filter(had_chemo == "Yes") %>% 
#   select(Chemotherapy = blood_bf_chemo_sample1) %>% 
#   tbl_summary()
# tbl2 <- breast_dna2 %>% #filter(had_hormone == "Yes") %>% 
#   select(Hormonetherapy = blood_bf_hormone_sample1) %>% 
#   tbl_summary()
# tbl3 <- breast_dna2 %>% #filter(had_immuno == "Yes") %>% 
#   select(Immnunotherapy = blood_bf_immuno_sample1) %>% 
#   tbl_summary()
# tbl4 <- breast_dna2 %>% #filter(had_rad == "Yes") %>% 
#   select(Radiotherapy = blood_bf_rad_sample1) %>% 
#   tbl_summary()
# tbl5 <- breast_dna2 %>% #filter(had_chemo_rad == "Yes") %>% 
#   select(`Chemotherapy+Radiation` = blood_bf_chemo_rad_sample1) %>% 
#   tbl_summary()
# tbl6 <- breast_dna2 %>% #filter(had_treatment == "Yes") %>% 
#   select(`All Treatment` = blood_bf_treatment_sample1) %>% 
#   tbl_summary()
# 
# tbl_s1 <- tbl_stack(list(tbl1, tbl2, tbl3, tbl4, tbl5, tbl6))
# 
# tbl1 <- breast_dna2 %>% filter(had_chemo == "Yes") %>% 
#   select(Chemotherapy = seq_sample_chemo) %>% 
#   tbl_summary()
# tbl2 <- breast_dna2 %>% filter(had_hormone == "Yes") %>% 
#   select(Hormonetherapy = seq_sample_hormone) %>% 
#   tbl_summary()
# tbl3 <- breast_dna2 %>% filter(had_immuno == "Yes") %>% 
#   select(Immnunotherapy = seq_sample_immuno) %>% 
#   tbl_summary()
# tbl4 <- breast_dna2 %>% filter(had_rad == "Yes") %>% 
#   select(Radiotherapy = seq_sample_rad) %>% 
#   tbl_summary()
# tbl5 <- breast_dna2 %>% filter(had_chemo_rad == "Yes") %>% 
#   select(`Chemotherapy+Radiation` = seq_sample_chemo) %>% 
#   tbl_summary()
# tbl6 <- breast_dna2 %>% filter(had_treatment == "Yes") %>% 
#   select(`All Treatment` = seq_sample_treatment) %>% 
#   tbl_summary()
# 
# tbl_s2 <- tbl_stack(list(tbl1, tbl2, tbl3, tbl4, tbl5, tbl6))
# 
# tbl_merge(list(tbl_s1, tbl_s2), tab_spanner = c("**Blood before**", "**Sequential blood after**")) %>% as_gt() %>%  
#   gt::tab_source_note(gt::md("**Unknown represent patients not treated for the corresponding category**"))


# Table if take samples -/+30 days treatment and sequential after treatment but with more tham 30 days after first blood????

















# breast_dna2 <- 
#   dcast(setDT(breast_dna1),
#         deidentified_patient_id+sample_family_id_sf+sample_id+specimen_collection_date+blood_bf_treatment+blood_bf_30_days_treatment ~ 
#           rowid(deidentified_patient_id),
#         value.var = c(
#                       "treatment_bf_blood", "treatment_bf_30days_blood", 
#                       "treatment_after_blood", "treatment_after_30days_blood")) %>% 
#   unite("treatment_bf_blood", starts_with("treatment_bf_blood"), sep = "; ", remove = TRUE, na.rm = TRUE) %>%
#   unite("treatment_bf_30days_blood", starts_with("treatment_bf_30days_blood"), sep = "; ", remove = TRUE, na.rm = TRUE) %>%
#   unite("treatment_after_blood", starts_with("treatment_after_blood"), sep = "; ", remove = TRUE, na.rm = TRUE) %>%
#   unite("treatment_after_30days_blood", starts_with("treatment_after_30days_blood"), sep = "; ", remove = TRUE, na.rm = TRUE) %>% 
#   mutate(across(where(is.character), ~ na_if(., ""))) %>% 
#   group_by(deidentified_patient_id) %>% 
#   mutate(sample_count = n()) %>% 
#   ungroup() %>% 
#   mutate(blood_bf_chemo = case_when(
#     treatment_bf_blood
#   ))

# Number of sample per patient
# breast_dna2 %>% distinct(deidentified_patient_id, .keep_all = TRUE) %>% 
#   select(sample_count) %>% 
#   tbl_summary(sort = list(everything() ~ "frequency"))
# 
# 
# breast_dna2 %>% 
#   ggplot(aes(x= blood_bf_treatment))+
#   geom_bar()+
#   coord_flip()
# 
# breast_dna3 <-
#   dcast(setDT(breast_dna2), deidentified_patient_id ~ rowid(deidentified_patient_id),
#         value.var = c("sample_family_id_sf", 
#                       "sample_id",
#                       "specimen_collection_date", 
#                       "blood_bf_treatment", "blood_bf_30_days_treatment", 
#                       "treatment_bf_blood", "treatment_bf_30days_blood",
#                       "treatment_after_blood", "treatment_after_30days_blood"),
#         sep = "_sample")


# breast_dna3 <- breast_dna2 %>% 
#   group_by(deidentified_patient_id) %>% 
#   summarise_at(vars(sample_family_id_sf, specimen_collection_date, blood_bf_treatment, blood_bf_30days_treatment, blood_after_treatment, blood_after_30days_treatment), paste, collapse = "/ ") 


# purrr::keep(~!all(is.na(.))) %>%
# select(deidentified_patient_id, ends_with("_1")) %>%
# `colnames<-`(str_remove(colnames(.), "_1"))


breast_info <- 
  breast_info %>% ##################################################### Fix bug
  left_join(., 
            Demographic %>% 
              select(deidentified_patient_id,
                     "date_of_birth"), 
            by = "date_of_birth")


# Second diagnosis
Second_dx <- Second_dx %>% 
  mutate(deidentified_patient_id = str_to_lower(deidentified_patient_id)) %>% 
  rename(cancer_site = primary_site)
# write_rds(Second_dx, "Second_dx.rds")

Second_dx %>% 
  select(cancer_site) %>% 
  tbl_summary(#sort = list(everything() ~ "frequency")
    )

Second_dx <- Second_dx %>% 
  # Create a variable for the first breast diagnose for each patient
  arrange(deidentified_patient_id, date_of_diagnosis) %>% 
  mutate(breast_cancer_dx = case_when(
    str_detect(primary_site, "BREAST")                  ~ date_of_diagnosis
  )) %>% 
  select(deidentified_patient_id, breast_cancer_dx, everything()) %>% 
  arrange(deidentified_patient_id, breast_cancer_dx) %>% 
  group_by(deidentified_patient_id) %>% 
  mutate(breast_cancer_dx = first(breast_cancer_dx)) %>% 
  ungroup() %>% 
  # Remove patients who never had breast cancer
  filter(!is.na(breast_cancer_dx)) %>% 
  # was the primary or secondary cancer
  mutate(relative_cancer = case_when(
    str_detect(primary_site, "BREAST") &
      breast_cancer_dx == date_of_diagnosis             ~ "breast",
    str_detect(primary_site, "BREAST")                  ~ "recidive",
    breast_cancer_dx < date_of_diagnosis                ~ "post-breast cancer",
    breast_cancer_dx > date_of_diagnosis                ~ "pre-breast cancer",
    breast_cancer_dx == date_of_diagnosis               ~ "diagnosed at the same time"
  ))

Second_dx %>% 
  filter(relative_cancer != "breast") %>% 
  select(cancer_site, relative_cancer) %>% 
  tbl_summary(by = relative_cancer)





################################################################################# III ### Merge data
Global_data <- full_join(Demographic, breast_dna, by = "deidentified_patient_id") %>% 
  full_join(., breast_info, by = "deidentified_patient_id") %>% 
  full_join(., Chemot, by = "deidentified_patient_id") %>% 
  full_join(., Hormonet, by = "deidentified_patient_id") %>% 
  full_join(., Immnunot, by = "deidentified_patient_id") %>% 
  full_join(., Radiot, by = "deidentified_patient_id")

write_rds(Global_data, "Global_data.rds")

blood_patients <- Global_data %>% 
  # filter to patients who have blood samples
  filter(!is.na(specimen_collection_date))

# End cleaning
