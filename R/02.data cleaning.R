################################################################################# II ### Data cleaning
# Demographic
Demographic <- Demographic %>%
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  # mutate(across(contains("date"), ~ as.Date(as.numeric(.), 
  #                                           origin = "1899-12-30")))
  mutate_at(c("mrn", "deidentified_patient_id"), ~str_to_lower(.))


# Chemotherapy
Chemot <- Chemot %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # Remove NA in both drug name and date
  filter_at(vars(chemotherapy_drug, chemotherapy_start_date,
                 chemotherapy_end_date), any_vars(!is.na(.))) %>% 
  # filter(!(is.na(chemotherapy_drug) & is.na(chemotherapy_start_date) & is.na(chemotherapy_end_date)))
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # Remove no chemo given in chemotherapy_drug
  filter(chemotherapy_drug != "no chemo given" | is.na(chemotherapy_drug))
  # # Remove no chemotherapy when 18.. or 2300 dates only but keep if real date
  # filter(!(str_detect(chemotherapy_start_date, "12:00:00 am") & # remove 271
  #            chemotherapy_completion_status_first_course == "no chemotherapy")) %>% 


# Check patient CONTRAINDICATED, NONE, NOT PLANNED, REFUSED
check <- Chemot %>% 
  filter(str_detect(chemotherapy_type, "contraindicated|none, not planned|refused|dc only")) %>% 
  arrange(chemotherapy_type)
write_csv(check, paste0(path, "/sanity check/chemo contraindicated|none, not planned|refused patients.csv"))
# write_csv(check %>% filter(str_detect(mrn, "271414|621813|772767|309843|748868")), paste0(path, "/sanity check/chemo contraindicated.csv"))
# write_csv(check %>% filter(str_detect(mrn, "571063|543238|1045598|729266|597917")), paste0(path, "/sanity check/chemo none, not planned.csv"))
# write_csv(check %>% filter(str_detect(mrn, "290870|307280|1034774|152952|176871")), paste0(path, "/sanity check/chemo refused.csv"))


Chemot1 <- Chemot %>% 
  # remove the chemo_type CONTRAINDICATED, PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(chemotherapy_type == "chemo nos" |
           chemotherapy_type == "multi-agent chemo" |
           chemotherapy_type == "none, not planned" | # clean more
           chemotherapy_type == "recommended,unkn if given" | # clean more
           chemotherapy_type == "single-agent chemo" |
           chemotherapy_type == "unknown; dc only"
  ) %>% 
  # Make it easier to not have na for future filtering, rescue the ones with a drug name
  mutate(chemotherapy_completion_status_first_course = case_when(
    !is.na(chemotherapy_drug) |
      !is.na(chemotherapy_start_date)                        ~ coalesce(chemotherapy_completion_status_first_course, "chemo given")
  ))

check <- Chemot1 %>% 
  filter(str_detect(chemotherapy_start_date, "12:00:00 am"))
write_csv(check, paste0(path, "/sanity check/chemo patients with 12:00:00 am dates.csv"))


Chemot1 <- Chemot1 %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(across(ends_with("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>% 
  # # To help after removing the unknown date
  # # mutate(drugs_unk = case_when(
  # #   str_detect(chemotherapy_start_date, "12:00:00 am") &
  # #     is.na(chemotherapy_drug)                              ~ "unk drug, 1800 date",
  # #   !is.na(chemotherapy_start_date) &
  # #     is.na(chemotherapy_drug)                              ~ "unk drug",
  # #   is.na(chemotherapy_start_date) &
  # #     is.na(chemotherapy_drug)                              ~ "unk drug, NA date"
  # # )) %>% 
  # mutate(data_unk = case_when(
  #   str_detect(chemotherapy_start_date, "12:00:00 am") &
  #     !is.na(chemotherapy_drug)                             ~ "known drug, 1800 date",
  #   !is.na(chemotherapy_start_date) &
  #     !is.na(chemotherapy_drug)                             ~ "known drug",
  #   is.na(chemotherapy_start_date) &
  #     !is.na(chemotherapy_drug)                             ~ "unk drug, NA date",
  #   str_detect(chemotherapy_start_date, "12:00:00 am") &
  #     is.na(chemotherapy_drug)                              ~ "unk drug, 1800 date",
  #   !is.na(chemotherapy_start_date) &
  #     is.na(chemotherapy_drug)                              ~ "unk drug, known date",
  # )) %>% 
  # # Dates
  # # remove the data 18.., 23..
  # mutate(chemotherapy_start_date = case_when(
  #   str_detect(chemotherapy_start_date, "12:00:00 am")       ~ NA_character_,
  #   TRUE                                                     ~ chemotherapy_start_date
  # ), 
  # chemotherapy_start_date = as.Date(as.numeric(chemotherapy_start_date), 
  #                                   origin = "1899-12-30")
  # ) %>% 
  # mutate(chemotherapy_end_date = case_when(
  #   str_detect(chemotherapy_end_date, "12:00:00 am")         ~ NA_character_,
  #   TRUE                                                     ~ chemotherapy_end_date
  # ), 
  # chemotherapy_end_date = as.Date(as.numeric(chemotherapy_end_date), 
  #                                 origin = "1899-12-30")
  # ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  # Fix the 2300 dates
mutate(chemotherapy_start_date = case_when(
  str_detect(chemotherapy_start_date, "2300")                   ~ NA_Date_,
  TRUE                                                          ~ chemotherapy_start_date
)) %>% 
  # mutate(treatment_unk = case_when(
  #   str_detect(chemotherapy_start_date, "2300|2301")          ~ "Unknown when and if given 2300",
  #   TRUE                                                      ~ NA_character_
  # )) %>% 
  mutate(chemotherapy_start_date =
           coalesce(chemotherapy_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%
  # mutate(chemotherapy_drug = coalesce(treatment_unk, chemotherapy_drug, data_unk)) %>% 
  # mutate(chemotherapy_start_date = case_when(
  #   str_detect(chemotherapy_start_date, "2300|2301")          ~ as.Date("1700-01-01", origin = "1899-12-30"),
  #   TRUE                                                      ~ chemotherapy_start_date
  # )) %>% 
  
  # Start pivot
  # select(mrn, chemotherapy_drug, chemotherapy_start_date, chemotherapy_end_date) %>% 
  distinct() 

Chemot <- Chemot1 %>%
  arrange(mrn, chemotherapy_start_date) %>% 
  # group_by(mrn, chemotherapy_start_date, chemotherapy_end_date) %>%
  group_by(mrn, chemotherapy_start_date) %>%
  # pivot_wider(id_cols = c(mrn, chemotherapy_start_date), values_from = chemotherapy_drug)
  summarise_at(vars(chemotherapy_drug, chemotherapy_end_date), str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(chemotherapy_end_date, paste("chemotherapy_end_date", 10:1, sep = ""), 
           sep = ";", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  purrr::keep(~!all(is.na(.)))

# pivot wider, use dcast bc better to keep date class
chemot <- dcast(setDT(Chemot), mrn ~ rowid(mrn),
                value.var = c(
                  "chemotherapy_drug",
                  "chemotherapy_start_date",
                  "chemotherapy_end_date1"
                ))

Chemot <- Chemot %>%
  select(mrn, treatment_start_date = chemotherapy_start_date,
         treatment_end_date = chemotherapy_end_date1, treatment = chemotherapy_drug) %>%
  mutate(treatment_type = "chemo")



# Hormonet
Hormonet <- Hormonet %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # remove rows with no drugs info
  filter_at(vars(hormone_therapy_drug, hormone_therapy_start_date,
                 hormone_therapy_end_date), any_vars(!is.na(.))) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # Remove no hormone given in drug
  filter(hormone_therapy_drug != "no hormone given" | is.na(hormone_therapy_drug))

# Check patient CONTRAINDICATED, NONE, NOT PLANNED, REFUSED
check <- Hormonet %>% 
  filter(str_detect(hormone_therapy_type, "contraindicated|none, not planned|refused|dc only")) %>% 
  arrange(hormone_therapy_type)
write_csv(check, paste0(path, "/sanity check/hormone contraindicated|none, not planned|refused patients.csv"))



Hormonet1 <- Hormonet %>% 
  # remove the chemo_type CONTRAINDICATED, PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(hormone_therapy_type == "hormone administered" |
           hormone_therapy_type == "none, not planned" |
           hormone_therapy_type == "recommended, unk if given" | # clean more
           hormone_therapy_type == "unknown; dc only"  # clean more
  )

check <- Hormonet1 %>% 
  filter(str_detect(hormone_therapy_start_date, "12:00:00 am"))
write_csv(check, paste0(path, "/sanity check/hormone patients with 12:00:00 am dates.csv"))

Hormonet1 <- Hormonet1 %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(across(ends_with("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>% 
  # mutate(data_unk = case_when(
  #   str_detect(hormone_therapy_start_date, "12:00:00 am") &
  #     !is.na(hormone_therapy_drug)                             ~ "known drug, 1800 date",
  #   !is.na(hormone_therapy_start_date) &
  #     !is.na(hormone_therapy_drug)                             ~ "known drug",
  #   is.na(hormone_therapy_start_date) &
  #     !is.na(hormone_therapy_drug)                             ~ "unk drug, NA date",
  #   str_detect(hormone_therapy_start_date, "12:00:00 am") &
  #     is.na(hormone_therapy_drug)                              ~ "unk drug, 1800 date",
  #   !is.na(hormone_therapy_start_date) &
  #     is.na(hormone_therapy_drug)                              ~ "unk drug, known date",
  # )) %>% 
  # filter(!is.na(data_unk)) %>% 
  # # Dates
  # # remove the data 18.., 23..
  # mutate(hormone_therapy_start_date = case_when(
  #   str_detect(hormone_therapy_start_date, "12:00:00 am")       ~ NA_character_,
  #   TRUE                                                     ~ hormone_therapy_start_date
  # ), 
  # hormone_therapy_start_date = as.Date(as.numeric(hormone_therapy_start_date), 
  #                                      origin = "1899-12-30")
  # ) %>% 
  # mutate(hormone_therapy_end_date = case_when(
  #   str_detect(hormone_therapy_end_date, "12:00:00 am")         ~ NA_character_,
  #   TRUE                                                     ~ hormone_therapy_end_date
  # ), 
  # hormone_therapy_end_date = as.Date(as.numeric(hormone_therapy_end_date), 
  #                                    origin = "1899-12-30")
  # ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  # Fix the 2300 dates
mutate(hormone_therapy_start_date = case_when(
  str_detect(hormone_therapy_start_date, "2300")                ~ NA_Date_,
  TRUE                                                          ~ hormone_therapy_start_date
)) %>% 
  # mutate(treatment_unk = case_when(
  #   str_detect(hormone_therapy_start_date, "2300|2301")          ~ "Unknown when and if given 2300",
  #   TRUE                                                      ~ NA_character_
  # )) %>% 
  mutate(hormone_therapy_start_date =
           coalesce(hormone_therapy_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # mutate(hormone_therapy_drug = coalesce(treatment_unk, hormone_therapy_drug, data_unk)) %>% 
  # mutate(hormone_therapy_start_date = case_when(
  #   str_detect(hormone_therapy_start_date, "2300|2301")          ~ as.Date("1700-01-01", origin = "1899-12-30"),
  #   TRUE                                                      ~ hormone_therapy_start_date
  # )) %>% 
  # 
  # # Start pivot
  # select(mrn, hormone_therapy_drug, hormone_therapy_start_date, hormone_therapy_end_date) %>% 
  distinct()

Hormonet <- Hormonet1 %>%
  arrange(mrn, hormone_therapy_start_date) %>% 
  # group_by(mrn, chemotherapy_start_date, chemotherapy_end_date) %>%
  group_by(mrn, hormone_therapy_start_date) %>%
  # pivot_wider(id_cols = c(mrn, hormone_therapy_start_date), values_from = chemotherapy_drug)
  summarise_at(vars(hormone_therapy_drug, hormone_therapy_end_date), str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(hormone_therapy_end_date, paste("hormone_therapy_end_date", 10:1, sep = ""), 
           sep = ";", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  purrr::keep(~!all(is.na(.)))


hormonet <- dcast(setDT(Hormonet), mrn ~ rowid(mrn),
                  value.var = c(
                    "hormone_therapy_drug",
                    "hormone_therapy_start_date",
                    "hormone_therapy_end_date1"
                  ))

Hormonet <- Hormonet %>%
  select(mrn, treatment_start_date = hormone_therapy_start_date,
         treatment_end_date = hormone_therapy_end_date1, treatment = hormone_therapy_drug) %>%
  mutate(treatment_type = "hormone")



# Immnunotherapy
Immnunot <- Immnunot %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # remove rows with no drugs info
  filter_at(vars(immunotherapy_drug, immunotherapy_start_date,
                 immunotherapy_end_date), any_vars(!is.na(.))) %>% 
  mutate(across(where(is.character), ~str_to_lower(.)))

# Check patient CONTRAINDICATED, NONE, NOT PLANNED, REFUSED
check <- Immnunot %>% 
  filter(str_detect(immunotherapy_type, "contraindicated|none, not planned|refused|dc only")) %>% 
  arrange(immunotherapy_type)
write_csv(check, paste0(path, "/sanity check/Immuno contraindicated|none, not planned|refused patients.csv"))

Immnunot1 <- Immnunot %>% 
  # remove the chemo_type CONTRAINDICATED, PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(immunotherapy_type == "immuno administered" |
           immunotherapy_type == "none, not planned" |
           immunotherapy_type == "recommended, unk if given"
  )

check <- Immnunot1 %>% 
  filter(str_detect(immunotherapy_start_date, "12:00:00 am"))
write_csv(check, paste0(path, "/sanity check/immuno patients with 12:00:00 am dates.csv"))
  
  
Immnunot1 <- Immnunot1 %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(across(ends_with("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>% 
  # mutate(immunotherapy_drug = case_when(
  #   is.na(immunotherapy_start_date) &
  #     !is.na(immunotherapy_drug)                             ~ "unk drug, NA date",
  #   str_detect(immunotherapy_start_date, "12:00:00 am") &
  #     is.na(immunotherapy_drug)                              ~ "unk drug, 1800 date",
  #   !is.na(immunotherapy_start_date) &
  #     is.na(immunotherapy_drug)                              ~ "unk drug, known date",
  # )) %>% 
  # filter(!is.na(immunotherapy_drug)) %>% 
  # # Dates
  # # remove the data 18.., 23..
  # mutate(immunotherapy_start_date = case_when(
  #   str_detect(immunotherapy_start_date, "12:00:00 am")       ~ NA_character_,
  #   TRUE                                                     ~ immunotherapy_start_date
  # ), 
  # immunotherapy_start_date = as.Date(as.numeric(immunotherapy_start_date), 
  #                                    origin = "1899-12-30")
  # ) %>% 
  # mutate(immunotherapy_end_date = case_when(
  #   str_detect(immunotherapy_end_date, "12:00:00 am")         ~ NA_character_,
  #   TRUE                                                     ~ immunotherapy_end_date
  # ), 
  # immunotherapy_end_date = as.Date(as.numeric(immunotherapy_end_date),
  #                                  origin = "1899-12-30")
  # ) %>% 
  # No 2300 date
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(immunotherapy_start_date =
           coalesce(immunotherapy_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # 
  # # Start pivot
  # select(mrn, immunotherapy_drug, immunotherapy_start_date, immunotherapy_end_date) %>% 
  distinct()

Immnunot <- Immnunot1 %>%
  arrange(mrn, immunotherapy_start_date) %>% 
  group_by(mrn, immunotherapy_start_date) %>%
  summarise_at(vars(immunotherapy_drug, immunotherapy_end_date), str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(immunotherapy_end_date, paste("immunotherapy_end_date", 10:1, sep = ""), 
           sep = ";", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  purrr::keep(~!all(is.na(.)))


immnunot <- dcast(setDT(Immnunot), mrn ~ rowid(mrn),
                  value.var = c(
                    "immunotherapy_drug",
                    "immunotherapy_start_date",
                    "immunotherapy_end_date1"
                  ))

Immnunot <- Immnunot %>%
  select(mrn, treatment_start_date = immunotherapy_start_date,
         treatment_end_date = immunotherapy_end_date1, treatment = immunotherapy_drug) %>%
  mutate(treatment_type = "immunoT")



# Radiotion
Radiot <- Radiot %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # remove rows with no rad info
  filter_at(vars(boost_dose_c_gy, radiation_start_date,
                 radiation_end_date), any_vars(!is.na(.))) %>% 
  mutate(across(where(is.character), ~str_to_lower(.)))

# Check patient CONTRAINDICATED, NONE, NOT PLANNED, REFUSED
check <- Radiot %>% 
  filter(str_detect(reason_for_no_radiation, "contraindicated|autopsy|refused|dco")) %>% 
  arrange(reason_for_no_radiation)
write_csv(check, paste0(path, "/sanity check/Rad contraindicated|autopsy|refused|dco patients.csv"))

Radiot1 <- Radiot %>% 
  # remove the chemo_type CONTRAINDICATED, PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(reason_for_no_radiation == "rad therapy performed" |
           # reason_for_no_radiation == "not recommended/autopsy" |
           reason_for_no_radiation == "recommended, unk if given" |
           reason_for_no_radiation == "unknown/dco"
  )

check <- Radiot1 %>% 
  filter(str_detect(radiation_start_date, "12:00:00 am"))
write_csv(check, paste0(path, "/sanity check/Rad patients with 12:00:00 am dates.csv"))


Radiot1 <- Radiot1 %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(across(ends_with("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>% 
  # mutate(data_unk = case_when(
  #   str_detect(radiation_start_date, "12:00:00 am") &
  #     !is.na(boost_dose_c_gy)                             ~ "known dose, 1800 date",
  #   !is.na(radiation_start_date) &
  #     !is.na(boost_dose_c_gy)                             ~ "known dose",
  #   is.na(radiation_start_date) &
  #     !is.na(boost_dose_c_gy)                             ~ "unk dose, NA date",
  #   str_detect(radiation_start_date, "12:00:00 am") &
  #     is.na(boost_dose_c_gy)                              ~ "unk dose, 1800 date",
  #   !is.na(radiation_start_date) &
  #     is.na(boost_dose_c_gy)                              ~ "unk dose, known date",
  # )) %>% 
  # filter(!is.na(data_unk)) %>% 
  # # Dates
  # # remove the data 18.., 23..
  # mutate(radiation_start_date = case_when(
  #   str_detect(radiation_start_date, "12:00:00 am")       ~ NA_character_,
  #   TRUE                                                     ~ radiation_start_date
  # ), 
  # radiation_start_date = as.Date(as.numeric(radiation_start_date), 
  #                                origin = "1899-12-30")
  # ) %>% 
  # mutate(radiation_end_date = case_when(
  #   str_detect(radiation_end_date, "12:00:00 am")         ~ NA_character_,
  #   TRUE                                                     ~ radiation_end_date
  # ), 
  # radiation_end_date = as.Date(as.numeric(radiation_end_date), ################################################ Fix bug create a 2300 date, lubridate::origin
  #                              origin = "1899-12-30")
  # ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  # Fix the 2300 dates
mutate(radiation_start_date = case_when(
  str_detect(radiation_start_date, "2300")                      ~ NA_Date_,
  TRUE                                                          ~ radiation_start_date
)) %>% 
  # mutate(treatment_unk = case_when(
  #   str_detect(radiation_start_date, "2300|2301")             ~ "Unknown when and if given 2300",
  #   TRUE                                                      ~ NA_character_
  # )) %>%
  mutate(radiation_start_date =
           coalesce(radiation_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # mutate(boost_dose_c_gy = coalesce(treatment_unk, as.character(boost_dose_c_gy), data_unk)) %>% 
  # mutate(radiation_start_date = case_when(
  #   str_detect(radiation_start_date, "2300|2301")             ~ as.Date("1700-01-01", origin = "1899-12-30"),
  #   TRUE                                                      ~ radiation_start_date
  # )) %>%
  # 
  # # Start pivot
  # select(mrn, boost_dose_c_gy, radiation_start_date, radiation_end_date) %>% 
  distinct()

Radiot <- Radiot1 %>%
  arrange(mrn, radiation_start_date) %>% 
  group_by(mrn, radiation_start_date) %>%
  summarise_at(vars(boost_dose_c_gy, radiation_end_date), str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(radiation_end_date, paste("radiation_end_date", 10:1, sep = ""), 
           sep = ";", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  purrr::keep(~!all(is.na(.)))



radiot <- dcast(setDT(Radiot), mrn ~ rowid(mrn),
                value.var = c(
                  "boost_dose_c_gy",
                  "radiation_start_date",
                  "radiation_end_date1"
                )
)

Radiot <- Radiot %>%
  select(mrn, treatment_start_date = radiation_start_date,
         treatment_end_date = radiation_end_date1, treatment = boost_dose_c_gy) %>%
  mutate(treatment_type = "radioT")




# Combine
treatment <- bind_rows(Chemot, Hormonet, Immnunot, Radiot) %>% 
  arrange(mrn, treatment_start_date) %>% 
  group_by(mrn, treatment_type) %>% 
  mutate(treatment_line = row_number(mrn)) %>% 
  unite(treatment_line, c(treatment_type, treatment_line), sep = "_", remove = FALSE)

Treatment <- full_join(chemot, hormonet, by = "mrn") %>% 
  full_join(., immnunot, by = "mrn") %>% 
  full_join(., radiot, by = "mrn") %>% 
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
  mutate(had_all_treatment = case_when(
    !is.na(chemotherapy_start_date_1) &
      !is.na(hormone_therapy_start_date_1) &
      !is.na(immunotherapy_start_date_1) &
      !is.na(radiation_start_date_1)             ~ "Yes"
  ))

write_rds(Treatment, "Treatment.rds")



# DNA, clean and add same sample/same date on the same row
table(breast_DNA$sample_type)

breast_dna <- breast_DNA %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  mutate_at(c("mrn"), ~str_to_lower(.)) %>% 
  mutate(across(contains("date"), ~ as.Date(as.numeric(.), 
                                            origin = "1899-12-30"))) %>% 
  filter(derived_tissue_type == "blood", str_detect(sample_type, "buffy coat|genomic dna|unprocessed liquid tissue|mnc less cd138+|mnc$|dna in prep")) %>% 
  mutate(deidentified_patient_id = str_to_lower(deidentified_patient_id)) %>% 
  select(mrn, deidentified_patient_id, sample_family_id_sf, sample_id,
         specimen_collection_date) %>%
  
  # mutate(specimen_collection_date = as.Date(specimen_collection_date, 
  #                                           format = "%Y-%M-%d", 
  #                                           origin = "1899-12-30")) %>% 
  arrange(mrn, specimen_collection_date) %>% 
  group_by(mrn, deidentified_patient_id, sample_family_id_sf, specimen_collection_date) %>% 
  summarise_at(vars(sample_id), str_c, collapse = "; ") %>%
  # separate(col = sample_id, paste("sample_id", 1:3, sep="_"), sep = "; ", extra = "drop", fill = "right")
  ungroup() #%>% 
  # left_join(., breast_DNA, 
  #           by = c("mrn", "deidentified_patient_id"))

write_rds(breast_dna, "breast_dna.rds")

# breast_dna1 <- breast_dna %>% left_join(., Treatment, by = "mrn") %>% 
#   mutate(blood_bf_chemo = case_when(
#     specimen_collection_date <= chemotherapy_start_date_1                ~ "Yes",
#     specimen_collection_date > chemotherapy_start_date_1                 ~ "No",
#     is.na(chemotherapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_hormone = case_when(
#     specimen_collection_date <= hormone_therapy_start_date_1                ~ "Yes",
#     specimen_collection_date > hormone_therapy_start_date_1                 ~ "No",
#     is.na(hormone_therapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_immuno = case_when(
#     specimen_collection_date <= immunotherapy_start_date_1                ~ "Yes",
#     specimen_collection_date > immunotherapy_start_date_1                 ~ "No",
#     is.na(immunotherapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_rad = case_when(
#     specimen_collection_date <= radiation_start_date_1                ~ "Yes",
#     specimen_collection_date > radiation_start_date_1                 ~ "No",
#     is.na(radiation_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_chemo_rad = case_when(
#     specimen_collection_date <= chemotherapy_start_date_1 &
#       specimen_collection_date <= radiation_start_date_1                ~ "Yes",
#     specimen_collection_date > chemotherapy_start_date_1 |
#       specimen_collection_date > radiation_start_date_1                 ~ "No",
#     is.na(chemotherapy_start_date_1) |
#       is.na(radiation_start_date_1)                                 ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_treatment = case_when(
#     specimen_collection_date <= chemotherapy_start_date_1 &
#       specimen_collection_date <= hormone_therapy_start_date_1 &
#       specimen_collection_date <= immunotherapy_start_date_1 &
#       specimen_collection_date <= radiation_start_date_1                ~ "Yes",
#     if_any(contains("start_date_1"), ~ . < specimen_collection_date)    ~ "No",
#     # if_all(contains("start_date_1"), ~ . < specimen_collection_date)    ~ "Nope",
#     is.na(chemotherapy_start_date_1) |
#       is.na(hormone_therapy_start_date_1) |
#       is.na(immunotherapy_start_date_1) |
#       is.na(radiation_start_date_1)                                     ~ "not administred",
#     TRUE                                                                ~ NA_character_
#   )) %>% 
#   mutate(blood_bf_30_days_chemo = case_when(
#     specimen_collection_date >= (chemotherapy_start_date_1 - days(30)) &
#       specimen_collection_date <= (chemotherapy_start_date_1 + days(30))              ~ "Yes",
#     specimen_collection_date > (chemotherapy_start_date_1 + days(30))                 ~ "No",
#     is.na(chemotherapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_hormone = case_when(
#     specimen_collection_date <= (hormone_therapy_start_date_1 + days(30))                ~ "Yes",
#     specimen_collection_date > (hormone_therapy_start_date_1 + days(30))                 ~ "No",
#     is.na(hormone_therapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_immuno = case_when(
#     specimen_collection_date <= (immunotherapy_start_date_1 + days(30))                ~ "Yes",
#     specimen_collection_date > (immunotherapy_start_date_1 + days(30))                 ~ "No",
#     is.na(immunotherapy_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_rad = case_when(
#     specimen_collection_date <= (radiation_start_date_1 + days(30))                ~ "Yes",
#     specimen_collection_date > (radiation_start_date_1 + days(30))                 ~ "No",
#     is.na(radiation_start_date_1)                                     ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_chemo_rad = case_when(
#     specimen_collection_date <= (chemotherapy_start_date_1 + days(30)) &
#       specimen_collection_date <= (radiation_start_date_1 + days(30))                ~ "Yes",
#     specimen_collection_date > (chemotherapy_start_date_1 + days(30)) &
#       specimen_collection_date > (radiation_start_date_1 + days(30))                 ~ "No",
#     is.na(chemotherapy_start_date_1) |
#       is.na(radiation_start_date_1)                                 ~ "not administred",
#     TRUE                                                            ~ NA_character_
#   )) %>%
#   mutate(blood_bf_30_days_treatment = case_when(
#     specimen_collection_date <= (chemotherapy_start_date_1 + days(30)) &
#       specimen_collection_date <= (hormone_therapy_start_date_1 + days(30)) &
#       specimen_collection_date <= (immunotherapy_start_date_1 + days(30)) &
#       specimen_collection_date <= (radiation_start_date_1 + days(30))                ~ "Yes",
#     if_any(contains("start_date_1"), ~ (.  + days(30)) < specimen_collection_date)    ~ "No",
#     # if_all(contains("start_date_1"), ~ . < specimen_collection_date)    ~ "Nope",
#     is.na(chemotherapy_start_date_1) |
#       is.na(hormone_therapy_start_date_1) |
#       is.na(immunotherapy_start_date_1) |
#       is.na(radiation_start_date_1)                                     ~ "not administred",
#     TRUE                                                                ~ NA_character_
#   )) %>%
#   mutate(across(contains("blood_bf_"), ~ factor(., levels = c("Yes", "No", "not administred")))) %>% 
#   
#   # ungroup() %>% 
#   # mutate(treatment_bf_blood = case_when(
#   #   specimen_collection_date <= treatment_start_date                ~ "Blood first",
#   #   specimen_collection_date <= treatment_start_date                ~ treatment_line,
#   #   TRUE                                                            ~ NA_character_
#   # )) %>% 
#   # mutate(treatment_bf_30days_blood = case_when(
#   #   specimen_collection_date <= (treatment_start_date + days(30))   ~ "Blood first",
#   #   specimen_collection_date <= (treatment_start_date + days(30))   ~ treatment_line,
#   #   TRUE                                                            ~ NA_character_
# # )) %>% 
# # 
# # mutate(treatment_after_blood = case_when(
# #   specimen_collection_date > treatment_start_date                 ~ treatment_line,
# #   TRUE                                                            ~ NA_character_
# # )) %>% 
# # mutate(treatment_after_30days_blood = case_when(
# #   specimen_collection_date > (treatment_start_date + days(30))   ~ treatment_line,
# #   TRUE                                                            ~ NA_character_
# # )) %>% 
# # distinct(mrn, sample_family_id_sf, sample_id,
# #          specimen_collection_date, )
# arrange(mrn, specimen_collection_date, blood_bf_treatment) %>% 
#   group_by(mrn) %>% 
#   # mutate(sample_lag = specimen_collection_date - lag(specimen_collection_date), 
#   #        sample_lag = str_remove(sample_lag, " days")
#   #        ) %>% 
#   # mutate(has_a_good_sample = case_when(
#   #   if_any(contains("blood_bf_"), ~ . == "Yes")    ~ "Yes"
#   # )) %>% 
#   # fill(has_a_good_sample, .direction = "down") %>% 
#   # mutate(has_a_good_seq_sample = case_when(
#   #   blood_bf_30_days_chemo_rad == "No" &
#   #     has_a_good_sample == "Yes" &
#   #     sample_lag > 30                             ~ "Yes"
# # )) %>% 
# 
# select(mrn, specimen_collection_date, #sample_lag, has_a_good_sample, has_a_good_seq_sample,
#        "blood_bf_chemo", "blood_bf_hormone", 
#        "blood_bf_immuno", "blood_bf_rad",
#        "blood_bf_chemo_rad",
#        "blood_bf_treatment", everything())
# 
# # breast_dna1 %>% distinct(mrn, .keep_all = TRUE) %>% 
# #   select(blood_bf_chemo, blood_bf_hormone, 
# #          blood_bf_immuno, blood_bf_rad,
# #          blood_bf_treatment, blood_bf_30_days_treatment) %>% 
# #   tbl_summary(sort = list(everything() ~ "frequency"))
# 
# breast_dna1 %>% filter(chemotherapy_start_date_1 == "1700-01-01") %>% nrow()
# breast_dna1 %>% filter(hormone_therapy_start_date_1 == "1700-01-01") %>% nrow()
# breast_dna1 %>% filter(immunotherapy_start_date_1 == "1700-01-01") %>% nrow()
# breast_dna1 %>% filter(radiation_start_date_1 == "1700-01-01") %>% nrow()
# breast_dna1 %>% filter(chemotherapy_start_date_1 == "1700-01-01" & 
#                          radiation_start_date_1 == "1700-01-01") %>% nrow()
# breast_dna1 %>% filter(chemotherapy_start_date_1 == "1700-01-01" & 
#                          hormone_therapy_start_date_1 == "1700-01-01" &
#                          immunotherapy_start_date_1 == "1700-01-01" &
#                          radiation_start_date_1 == "1700-01-01") %>% nrow()
# 
# 
# 
# # For samples before treatment, take only the closest
# 
# sample_before_chemo <- breast_dna1 %>% 
#   filter(blood_bf_chemo == "Yes") %>% 
#   # pivot_longer(cols = c(chemotherapy_start_date_1, hormone_therapy_start_date_1, 
#   #                       immunotherapy_start_date_1, radiation_start_date_1),
#   #              names_to = "treatmetn_type", values_to = "date_treatment") %>% 
#   mutate(interval_sample_chemo = 
#            abs(interval(start = specimen_collection_date, end = chemotherapy_start_date_1) /
#                  duration(n = 1, units = "days"))) %>% 
#   arrange(mrn, interval_sample_chemo) %>% 
#   distinct(mrn, .keep_all = TRUE)
# 
# library(ggforce)
# p <- qplot(x =interval_sample_chemo, data=subset(sample_before_chemo), fill=..count.., 
#            geom="histogram", 
#            binwidth = 100,
# ) 
# p + scale_fill_viridis_c(
#   alpha = 1,
#   begin = 0,
#   end = 1,
#   direction = 1,
#   option = "D",
#   values = NULL,
#   space = "Lab",
#   na.value = "grey50",
#   guide = "colourbar",
#   aesthetics = "fill"
# ) +
#   theme_minimal(base_size = 14) +
#   labs(x="Time from Blood Collection to Chemotherapy (in days)", 
#        y="Number of Patient",
#        caption = "Each bar represents 100 days") +
#   facet_zoom(ylim = c(0, 10), zoom.size = 1)
# 
# sample_before_hormone <- breast_dna1 %>% 
#   filter(blood_bf_hormone == "Yes") %>% 
#   mutate(interval_sample_hormone = 
#            abs(interval(start = specimen_collection_date, end = hormone_therapy_start_date_1) /
#                  duration(n = 1, units = "days"))) %>% 
#   arrange(mrn, interval_sample_hormone) %>% 
#   distinct(mrn, .keep_all = TRUE)
# 
# p <- qplot(x =interval_sample_hormone, data=subset(sample_before_hormone), fill=..count.., 
#            geom="histogram", 
#            binwidth = 100,
# ) 
# p + scale_fill_viridis_c(
#   alpha = 1,
#   begin = 0,
#   end = 1,
#   direction = 1,
#   option = "A",
#   values = NULL,
#   space = "Lab",
#   na.value = "grey50",
#   guide = "colourbar",
#   aesthetics = "fill"
# ) +
#   theme_minimal(base_size = 14) +
#   labs(x="Time from Blood Collection to Hormonetherapy (in days)", 
#        y="Number of Patient",
#        caption = "Each bar represents 100 days") +
#   facet_zoom(ylim = c(0, 10), zoom.size = 1)
# 
# sample_before_immuno <- breast_dna1 %>% 
#   filter(blood_bf_immuno == "Yes") %>% 
#   mutate(interval_sample_immuno = 
#            abs(interval(start = specimen_collection_date, end = immunotherapy_start_date_1) /
#                  duration(n = 1, units = "days"))) %>% 
#   arrange(mrn, interval_sample_immuno) %>% 
#   distinct(mrn, .keep_all = TRUE)
# 
# p <- qplot(x =interval_sample_immuno, data=subset(sample_before_immuno), fill=..count.., 
#            geom="histogram", 
#            binwidth = 100,
# ) 
# p + scale_fill_viridis_c(
#   alpha = 1,
#   begin = 0,
#   end = 1,
#   direction = 1,
#   option = "H",
#   values = NULL,
#   space = "Lab",
#   na.value = "grey50",
#   guide = "colourbar",
#   aesthetics = "fill"
# ) +
#   theme_minimal(base_size = 14) +
#   labs(x="Time from Blood Collection to Immunotherapy (in days)", 
#        y="Number of Patient",
#        caption = "Each bar represents 100 days") +
#   facet_zoom(ylim = c(0, 10), zoom.size = 1)
# 
# sample_before_rad <- breast_dna1 %>% 
#   filter(blood_bf_rad == "Yes") %>% 
#   mutate(interval_sample_rad = 
#            abs(interval(start = specimen_collection_date, end = radiation_start_date_1) /
#                  duration(n = 1, units = "days"))) %>% 
#   arrange(mrn, interval_sample_rad) %>% 
#   distinct(mrn, .keep_all = TRUE)
# 
# p <- qplot(x =interval_sample_rad, data=subset(sample_before_rad), fill=..count.., 
#            geom="histogram", 
#            binwidth = 100,
# ) 
# p + scale_fill_viridis_c(
#   alpha = 1,
#   begin = 0,
#   end = 1,
#   direction = 1,
#   option = "A",
#   values = NULL,
#   space = "Lab",
#   na.value = "grey50",
#   guide = "colourbar",
#   aesthetics = "fill"
# ) +
#   theme_minimal(base_size = 14) +
#   labs(x="Time from Blood Collection to Radiotherapy (in days)", 
#        y="Number of Patient",
#        caption = "Each bar represents 100 days") +
#   facet_zoom(ylim = c(0, 5), zoom.size = 1)
# 
# sample_before <- bind_rows(sample_before_chemo, sample_before_hormone, sample_before_immuno, sample_before_rad)
# 
# 
# 
# breast_dna2 <- breast_dna1 %>% 
#   mutate(had_good_sample_chemo = case_when(
#     blood_bf_chemo == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(mrn) %>% 
#   fill(had_good_sample_chemo, .direction = "updown") %>% 
#   mutate(seq_sample_chemo = case_when(
#     had_good_sample_chemo == "Yes" &
#       blood_bf_chemo == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_hormone = case_when(
#     blood_bf_hormone == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(mrn) %>% 
#   fill(had_good_sample_hormone, .direction = "updown") %>% 
#   mutate(seq_sample_hormone = case_when(
#     had_good_sample_hormone == "Yes" &
#       blood_bf_hormone == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_immuno = case_when(
#     blood_bf_immuno == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(mrn) %>% 
#   fill(had_good_sample_immuno, .direction = "updown") %>% 
#   mutate(seq_sample_immuno = case_when(
#     had_good_sample_immuno == "Yes" &
#       blood_bf_immuno == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_rad = case_when(
#     blood_bf_rad == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(mrn) %>% 
#   fill(had_good_sample_rad, .direction = "updown") %>% 
#   mutate(seq_sample_rad = case_when(
#     had_good_sample_rad == "Yes" &
#       blood_bf_rad == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_chemo_rad = case_when(
#     blood_bf_chemo_rad == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(mrn) %>% 
#   fill(had_good_sample_chemo_rad, .direction = "updown") %>% 
#   mutate(seq_sample_chemo_rad = case_when(
#     had_good_sample_chemo_rad == "Yes" &
#       blood_bf_chemo_rad == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   mutate(had_good_sample_treatment = case_when(
#     blood_bf_treatment == "Yes"              ~ "Yes"
#   )) %>% 
#   group_by(mrn) %>% 
#   fill(had_good_sample_treatment, .direction = "updown") %>% 
#   mutate(seq_sample_treatment = case_when(
#     had_good_sample_treatment == "Yes" &
#       blood_bf_treatment == "No" ~ "Yes",
#     TRUE ~ "No"
#   )) %>% 
#   ungroup()
# 
# 
# 
# 
# # left_join(., Treatment %>% 
# #             select(mrn, had_chemo, had_hormone, had_immuno, had_rad, had_chemo_rad, had_treatment),
# #           by = "mrn")
# 
# sample_after_chemo <- breast_dna2 %>%
#   filter(seq_sample_chemo == "Yes") %>%
#   mutate(interval_sample_chemo =
#            abs(interval(start = specimen_collection_date, end = chemotherapy_start_date_1) /
#                  duration(n = 1, units = "days"))) %>%
#   arrange(mrn, interval_sample_chemo) %>%
#   distinct(mrn, specimen_collection_date, .keep_all = TRUE) %>% 
#   group_by(mrn) %>% 
#   mutate(sequential_sample_count = factor(row_number(mrn))) %>% 
#   ungroup() %>% 
#   select(mrn, sequential_sample_count, interval_sample_chemo, everything())
# 
# sample_after_chemo %>% 
#   ggplot(aes(x =interval_sample_chemo, fill = sequential_sample_count))+
#   geom_histogram(binwidth = 100, alpha = 0.9, position = "stack") +
#   scale_fill_viridis(discrete=T)+
#   theme_minimal(base_size = 14) +
#   labs(x="Time from Blood Collection to Chemotherapy (in days)", 
#        y="Number of Patient",
#        caption = "Each bar represents 100 days")
# 
# sample_after_hormone <- breast_dna2 %>%
#   filter(seq_sample_hormone == "Yes") %>%
#   mutate(interval_sample_hormone =
#            abs(interval(start = specimen_collection_date, end = hormone_therapy_start_date_1) /
#                  duration(n = 1, units = "days"))) %>%
#   arrange(mrn, interval_sample_hormone) %>%
#   distinct(mrn, specimen_collection_date, .keep_all = TRUE) %>% 
#   group_by(mrn) %>% 
#   mutate(sequential_sample_count = factor(row_number(mrn))) %>% 
#   ungroup() %>% 
#   select(mrn, sequential_sample_count, interval_sample_hormone, everything())
# 
# sample_after_hormone %>% 
#   ggplot(aes(x =interval_sample_hormone, fill = sequential_sample_count))+
#   geom_histogram(binwidth = 100, alpha = 0.9, position = "stack") +
#   scale_fill_viridis(discrete=T)+
#   theme_minimal(base_size = 14) +
#   labs(x="Time from Blood Collection to Hormonetherapy (in days)", 
#        y="Number of Patient",
#        caption = "Each bar represents 100 days")
# 
# sample_after_immuno <- breast_dna2 %>%
#   filter(seq_sample_immuno == "Yes") %>%
#   mutate(interval_sample_immuno =
#            abs(interval(start = specimen_collection_date, end = immunotherapy_start_date_1) /
#                  duration(n = 1, units = "days"))) %>%
#   arrange(mrn, interval_sample_immuno) %>%
#   distinct(mrn, specimen_collection_date, .keep_all = TRUE) %>% 
#   group_by(mrn) %>% 
#   mutate(sequential_sample_count = factor(row_number(mrn))) %>% 
#   ungroup() %>% 
#   select(mrn, sequential_sample_count, interval_sample_immuno, everything())
# 
# sample_after_immuno %>% 
#   ggplot(aes(x =interval_sample_immuno, fill = sequential_sample_count))+
#   geom_histogram(binwidth = 100, alpha = 0.9, position = "stack") +
#   scale_fill_viridis(discrete=T)+
#   theme_minimal(base_size = 14) +
#   labs(x="Time from Blood Collection to Immunotherapy (in days)", 
#        y="Number of Patient",
#        caption = "Each bar represents 100 days")
# 
# sample_after_rad <- breast_dna2 %>%
#   filter(seq_sample_rad == "Yes") %>%
#   mutate(interval_sample_rad =
#            abs(interval(start = specimen_collection_date, end = radiation_start_date_1) /
#                  duration(n = 1, units = "days"))) %>%
#   arrange(mrn, interval_sample_rad) %>%
#   distinct(mrn, specimen_collection_date, .keep_all = TRUE) %>% 
#   group_by(mrn) %>% 
#   mutate(sequential_sample_count = factor(row_number(mrn))) %>% 
#   ungroup() %>% 
#   select(mrn, sequential_sample_count, interval_sample_rad, everything())
# 
# sample_after_rad %>% 
#   ggplot(aes(x =interval_sample_rad, fill = sequential_sample_count))+
#   geom_histogram(binwidth = 100, alpha = 0.9, position = "stack") +
#   scale_fill_viridis(discrete=T)+
#   theme_minimal(base_size = 14) +
#   labs(x="Time from Blood Collection to Radiotherapy (in days)", 
#        y="Number of Patient",
#        caption = "Each bar represents 100 days")
# 
# sample_after <- bind_rows(sample_after_chemo, sample_after_hormone, sample_after_immuno, sample_after_rad)













# breast_dna2 <-
#   dcast(setDT(breast_dna1), mrn ~ rowid(mrn),
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
#               select(mrn, had_chemo, had_hormone, had_immuno, had_rad, had_chemo_rad, had_treatment),
#             by = "mrn") %>% 
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
#         mrn+sample_family_id_sf+sample_id+specimen_collection_date+blood_bf_treatment+blood_bf_30_days_treatment ~ 
#           rowid(mrn),
#         value.var = c(
#                       "treatment_bf_blood", "treatment_bf_30days_blood", 
#                       "treatment_after_blood", "treatment_after_30days_blood")) %>% 
#   unite("treatment_bf_blood", starts_with("treatment_bf_blood"), sep = "; ", remove = TRUE, na.rm = TRUE) %>%
#   unite("treatment_bf_30days_blood", starts_with("treatment_bf_30days_blood"), sep = "; ", remove = TRUE, na.rm = TRUE) %>%
#   unite("treatment_after_blood", starts_with("treatment_after_blood"), sep = "; ", remove = TRUE, na.rm = TRUE) %>%
#   unite("treatment_after_30days_blood", starts_with("treatment_after_30days_blood"), sep = "; ", remove = TRUE, na.rm = TRUE) %>% 
#   mutate(across(where(is.character), ~ na_if(., ""))) %>% 
#   group_by(mrn) %>% 
#   mutate(sample_count = n()) %>% 
#   ungroup() %>% 
#   mutate(blood_bf_chemo = case_when(
#     treatment_bf_blood
#   ))

# Number of sample per patient
# breast_dna2 %>% distinct(mrn, .keep_all = TRUE) %>% 
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
#   dcast(setDT(breast_dna2), mrn ~ rowid(mrn),
#         value.var = c("sample_family_id_sf", 
#                       "sample_id",
#                       "specimen_collection_date", 
#                       "blood_bf_treatment", "blood_bf_30_days_treatment", 
#                       "treatment_bf_blood", "treatment_bf_30days_blood",
#                       "treatment_after_blood", "treatment_after_30days_blood"),
#         sep = "_sample")


# breast_dna3 <- breast_dna2 %>% 
#   group_by(mrn) %>% 
#   summarise_at(vars(sample_family_id_sf, specimen_collection_date, blood_bf_treatment, blood_bf_30days_treatment, blood_after_treatment, blood_after_30days_treatment), paste, collapse = "/ ") 


# purrr::keep(~!all(is.na(.))) %>%
# select(mrn, ends_with("_1")) %>%
# `colnames<-`(str_remove(colnames(.), "_1"))


breast_info <- 
  breast_info %>%
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  mutate_at(c("mrn"), ~str_to_lower(.)) %>% 
  filter(str_detect(primary_site, "breast")) %>% 
  distinct(mrn, date_of_diagnosis, .keep_all = TRUE) %>% 
  mutate(across(contains("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>%
  arrange(mrn, date_of_diagnosis)
  # left_join(., 
  #           Demographic %>% 
  #             select(mrn,
  #                    "date_of_birth"), 
  #           by = c("mrn", "date_of_birth"))
breast_info1 <- dcast(setDT(breast_info), mrn+date_of_birth ~ rowid(mrn),
                value.var = c(
                  "date_of_diagnosis",
                  "primary_site",
                  "histology",
                  "first_treatment_date",
                  "summary_of_rx_1st_course"
                )
)

# # Second diagnosis
# Second_dx <- Second_dx %>% 
#   mutate(mrn = str_to_lower(mrn)) %>% 
#   rename(cancer_site = primary_site)
# # write_rds(Second_dx, "Second_dx.rds")
# 
# Second_dx %>% 
#   select(cancer_site) %>% 
#   tbl_summary(#sort = list(everything() ~ "frequency")
#     )
# 
# Second_dx <- Second_dx %>% 
#   # Create a variable for the first breast diagnose for each patient
#   arrange(mrn, date_of_diagnosis) %>% 
#   mutate(breast_cancer_dx = case_when(
#     str_detect(primary_site, "BREAST")                  ~ date_of_diagnosis
#   )) %>% 
#   select(mrn, breast_cancer_dx, everything()) %>% 
#   arrange(mrn, breast_cancer_dx) %>% 
#   group_by(mrn) %>% 
#   mutate(breast_cancer_dx = first(breast_cancer_dx)) %>% 
#   ungroup() %>% 
#   # Remove patients who never had breast cancer
#   filter(!is.na(breast_cancer_dx)) %>% 
#   # was the primary or secondary cancer
#   mutate(relative_cancer = case_when(
#     str_detect(primary_site, "BREAST") &
#       breast_cancer_dx == date_of_diagnosis             ~ "breast",
#     str_detect(primary_site, "BREAST")                  ~ "recidive",
#     breast_cancer_dx < date_of_diagnosis                ~ "post-breast cancer",
#     breast_cancer_dx > date_of_diagnosis                ~ "pre-breast cancer",
#     breast_cancer_dx == date_of_diagnosis               ~ "diagnosed at the same time"
#   ))
# 
# Second_dx %>% 
#   filter(relative_cancer != "breast") %>% 
#   select(cancer_site, relative_cancer) %>% 
#   tbl_summary(by = relative_cancer)





################################################################################# III ### Merge data
Global_data <- full_join(Demographic, breast_dna, by = c("mrn", "deidentified_patient_id")) %>% 
  full_join(., breast_info1, by = c("mrn", "date_of_birth")) %>% 
  full_join(., chemot, by = "mrn") %>% 
  full_join(., hormonet, by = "mrn") %>% 
  full_join(., immnunot, by = "mrn") %>% 
  full_join(., radiot, by = "mrn")

write_rds(Global_data, "Global_data.rds")


# End cleaning
