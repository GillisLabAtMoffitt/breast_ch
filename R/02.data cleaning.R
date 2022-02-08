################################################################################# II ### Data cleaning

# Breast samples----
table(breast_DNA$sample_type)

breast_dna <- breast_DNA %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  mutate_at(c("mrn"), ~str_to_lower(.)) %>% 
  filter(collection_site_tissue_type == "Blood", 
         str_detect(sample_type, "Buffy|Genomic|Unprocessed|CD138|MNC$")) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  select(mrn, party_id, sample_family_id, sample_id,
         specimen_collection_date) %>%
  # add same sample/same date on the same row
  arrange(mrn, specimen_collection_date) %>% 
  group_by(mrn, party_id, sample_family_id, specimen_collection_date) %>% 
  summarise_at(vars(sample_id), str_c, collapse = "; ") %>%
  # separate(col = sample_id, paste("sample_id", 1:3, sep="_"), sep = "; ", extra = "drop", fill = "right")
  ungroup()

write_rds(breast_dna, "breast_dna.rds")


# Cancer Characteristics----
breast_info <- 
  breast_info %>%
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)
         ) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  filter(str_detect(primary_site_group, "Breast")) %>%
  distinct(mrn, date_of_diagnosis, .keep_all = TRUE) %>% 
  mutate_at("date_of_diagnosis", ~ as.Date(as.numeric(.), 
                                            origin = "1899-12-30")) %>%
  arrange(mrn, date_of_diagnosis) %>% 
  select(-c(group_name))

# breast_info1 <- dcast(setDT(breast_info), mrn+date_of_birth ~ rowid(mrn),
#                       value.var = c(
#                         "date_of_diagnosis",
#                         "primary_site",
#                         "class_of_case",
#                         "histology",
#                         "clinical_tnm_group_stage",
#                         "summary_of_rx_1st_course"
#                       )) %>% 
#   purrr::keep(~!all(is.na(.)))


# Demographic----
Demographic <- Demographic %>%
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  select(-c(group_name))


################################################################################# III ### Merge data
breast_patients <- breast_dna %>% 
  # Merge with Cancer Char for patients with samples available
  left_join(., breast_info, by = c("mrn", "party_id")) %>% 
  # Merge with Demographic for patients with samples available
  left_join(., Demographic, 
            c("mrn", "party_id", "date_of_birth"))

breast_patients_id <- paste0(breast_patients$mrn, collapse = "|")


################################################################################# IV ### Clean treatments
# Chemotherapy----
Chemot <- Chemot %>% 
  filter(!str_detect(chemotherapy_drug, "AG-013736|BRENTUXUMAB|DEPOCYT") |
         is.na(chemotherapy_drug)) %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # Limit to Breast cancer patients
  filter(str_detect(mrn, breast_patients_id)) %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(chemotherapy_end_date = as.Date(as.numeric(chemotherapy_end_date), 
                                             origin = "1899-12-30")) %>% 
  mutate(chemotherapy_start_date = as.Date(chemotherapy_start_date)) %>% 
  # Fix the 2300 dates
  mutate(chemotherapy_start_date = case_when(
    str_detect(chemotherapy_start_date, "2300")                   ~ NA_Date_,
    TRUE                                                          ~ chemotherapy_start_date
  )) %>% 
  mutate(chemotherapy_end_date = case_when(
    str_detect(chemotherapy_end_date, "2300")                     ~ NA_Date_,
    TRUE                                                          ~ chemotherapy_end_date
  )) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # Remove no chemo given in chemotherapy_drug
  filter(chemotherapy_drug != "no chemo given" | is.na(chemotherapy_drug)) %>% 
  # Remove NA in both drug name and date
  filter_at(vars(chemotherapy_drug, chemotherapy_start_date,
                 chemotherapy_end_date), any_vars(!is.na(.)))
  # filter(!(is.na(chemotherapy_drug) & is.na(chemotherapy_start_date) & is.na(chemotherapy_end_date)))
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
           chemotherapy_type == "contraindicated" |
           chemotherapy_type == "none, not planned" | # clean more
           chemotherapy_type == "recommended,unkn if given" |
           chemotherapy_type == "single-agent chemo" |
           chemotherapy_type == "unknown; dc only"
  ) %>% 
  # Make it easier to not have na for future filtering, rescue the ones with a drug name
  mutate(chemotherapy_completion_status_first_course = case_when(
    !is.na(chemotherapy_drug) |
      !is.na(chemotherapy_start_date)         
    ~ coalesce(chemotherapy_completion_status_first_course, "chemo given")
  ))

check <- Chemot1 %>%
  filter(str_detect(chemotherapy_completion_status_first_course, "no chem"))
write_csv(check, paste0(path, "/sanity check/chemo patients no chemo given.csv"))


Chemot1 <- Chemot1 %>% 
  # Clean contraindicated and none, not planned
  mutate(remove = case_when(
    chemotherapy_type == "contraindicated" &
      is.na(chemotherapy_drug)                                    ~ 1,
    chemotherapy_type == "none, not planned" &
      (is.na(chemotherapy_drug) & 
         is.na(chemotherapy_start_date))                          ~ 1,
    TRUE                                                          ~ 0
  )) %>% 
  filter(remove == 0) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(chemotherapy_start_date =
           coalesce(chemotherapy_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%

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

  # mutate(treatment_unk = case_when(
  #   str_detect(chemotherapy_start_date, "2300|2301")          ~ "Unknown when and if given 2300",
  #   TRUE                                                      ~ NA_character_
  # )) %>% 
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
  summarise_at(vars(chemotherapy_drug, chemotherapy_end_date), str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(chemotherapy_end_date, paste("chemotherapy_end_date", 10:1, sep = ""), 
           sep = "; ", remove = TRUE, 
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



# Hormonetherapy----
Hormonet <- Hormonet %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # Limit to Breast cancer patients
  filter(str_detect(mrn, breast_patients_id)) %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(across(ends_with("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>% 
  # Fix the 2300 dates
  mutate(hormone_therapy_start_date = case_when(
    str_detect(hormone_therapy_start_date, "2300")                ~ NA_Date_,
    TRUE                                                          ~ hormone_therapy_start_date
  )) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # Remove no hormone given in hormone_therapy_drug
  filter(hormone_therapy_drug != "no hormone given" | is.na(hormone_therapy_drug)) %>% 
  # Remove NA in both drug name and date
  filter_at(vars(hormone_therapy_drug, hormone_therapy_start_date,
                 hormone_therapy_end_date), any_vars(!is.na(.)))

# Check patient CONTRAINDICATED, NONE, NOT PLANNED, REFUSED
check <- Hormonet %>% 
  filter(str_detect(hormone_therapy_type, "contraindicated|none, not planned|refused|dc only")) %>% 
  arrange(hormone_therapy_type)
write_csv(check, paste0(path, "/sanity check/hormone contraindicated|none, not planned|refused patients.csv"))

Hormonet1 <- Hormonet %>% 
  # remove the chemo_type, PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(hormone_therapy_type == "hormone administered" |
           hormone_therapy_type == "none, not planned" | # clean more
           hormone_therapy_type == "contraindicated" |
           hormone_therapy_type == "recommended, unk if given" | 
           hormone_therapy_type == "unknown; dc only" 
  )

# check <- Hormonet1 %>% 
#   filter(str_detect(hormone_therapy_start_date, "12:00:00 am"))
# write_csv(check, paste0(path, "/sanity check/hormone patients with 12:00:00 am dates.csv"))

Hormonet1 <- Hormonet1 %>% 
  # Clean contraindicated and none, not planned
  mutate(remove = case_when(
    hormone_therapy_type == "contraindicated" &
      (is.na(hormone_therapy_drug) & 
         is.na(hormone_therapy_start_date))                       ~ 1,
    hormone_therapy_type == "none, not planned" &
      (is.na(hormone_therapy_drug) & 
         is.na(hormone_therapy_start_date))                       ~ 1,
    TRUE                                                          ~ 0
  )) %>% 
  filter(remove == 0) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(hormone_therapy_start_date =
           coalesce(hormone_therapy_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% 
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
  # mutate(treatment_unk = case_when(
  #   str_detect(hormone_therapy_start_date, "2300|2301")          ~ "Unknown when and if given 2300",
  #   TRUE                                                      ~ NA_character_
  # )) %>% 
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



# Immnunotherapy----
Immnunot <- Immnunot %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # Limit to Breast cancer patients
  filter(str_detect(mrn, breast_patients_id)) %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(across(ends_with("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>% 
  # No 2300 date
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # Remove no hormone given in immunotherapy_drug
  # filter(immunotherapy_drug != "no  given" | is.na(immunotherapy_drug)) %>% 
  # remove rows with no drugs info
  filter_at(vars(immunotherapy_drug, immunotherapy_start_date,
                 immunotherapy_end_date), any_vars(!is.na(.)))

# Check patient CONTRAINDICATED, NONE, NOT PLANNED, REFUSED
check <- Immnunot %>% 
  filter(str_detect(immunotherapy_type, "contraindicated|none, not planned|refused|dc only")) %>% 
  arrange(immunotherapy_type)
write_csv(check, paste0(path, "/sanity check/Immuno contraindicated|none, not planned|refused patients.csv"))

Immnunot1 <- Immnunot %>% 
  # remove the chemo_type CONTRAINDICATED, PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(immunotherapy_type == "immuno administered" |
           immunotherapy_type == "none, not planned" |
           immunotherapy_type == "contraindicated" |
           immunotherapy_type == "unknown; dc only" |
           immunotherapy_type == "recommended, unk if given"
  )

# check <- Immnunot1 %>% 
#   filter(str_detect(immunotherapy_start_date, "12:00:00 am"))
# write_csv(check, paste0(path, "/sanity check/immuno patients with 12:00:00 am dates.csv"))


Immnunot1 <- Immnunot1 %>% 
  # Clean contraindicated and none, not planned
  mutate(remove = case_when(
    immunotherapy_type == "contraindicated" &
      (is.na(immunotherapy_drug) & 
         is.na(immunotherapy_start_date))                        ~ 1,
    immunotherapy_type == "none, not planned" &
      (is.na(immunotherapy_drug) & 
         is.na(immunotherapy_start_date))                        ~ 1,
    TRUE                                                         ~ 0
  )) %>% 
  filter(remove == 0) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(immunotherapy_start_date =
           coalesce(immunotherapy_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%

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



# Radiotion----
Radiot <- Radiot %>% 
  mutate(mrn = as.character(mrn),
         mrn = coalesce(mrn, party_id)) %>% 
  # Limit to Breast cancer patients
  filter(str_detect(mrn, breast_patients_id)) %>% 
  # Transform number to date and 12:00:00 am character as NA
  mutate(across(ends_with("date"), ~ as.Date(as.numeric(.), 
                                             origin = "1899-12-30"))) %>% 
  # Fix the 2300 dates
  mutate(radiation_start_date = case_when(
    str_detect(radiation_start_date, "2300")                      ~ NA_Date_,
    TRUE                                                          ~ radiation_start_date
  )) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  # remove rows with no rad info
  filter(!is.na(boost_dose_c_gy) | !is.na(radiation_start_date),
         !is.na(radiation_end_date) | 
         radiation_location_of_rx != "no radiation therapy")

# Check patient CONTRAINDICATED, NONE, NOT PLANNED, REFUSED
check <- Radiot %>% 
  filter(str_detect(reason_for_no_radiation, "contraindicated|autopsy|refused|dco")) %>% 
  arrange(reason_for_no_radiation)
write_csv(check, paste0(path, "/sanity check/Rad contraindicated|autopsy|refused|dco patients.csv"))

Radiot1 <- Radiot %>% 
  # remove the PT DIED, RECOMMENDED, NOT GIVEN, "REFUSED"
  filter(reason_for_no_radiation == "rad therapy performed" |
           # reason_for_no_radiation == "not recommended/autopsy" |
           reason_for_no_radiation == "recommended, unk if given" |
           reason_for_no_radiation == "radiation contraindicated" |
           reason_for_no_radiation == "unknown/dco"
  )

# check <- Radiot1 %>% 
#   filter(str_detect(radiation_start_date, "12:00:00 am"))
# write_csv(check, paste0(path, "/sanity check/Rad patients with 12:00:00 am dates.csv"))


Radiot1 <- Radiot1 %>% 
  # Clean contraindicated and none, not planned
  mutate(remove = case_when(
    reason_for_no_radiation == "radiation contraindicated" &
      is.na(radiation_start_date)                                 ~ 1,
    TRUE                                                          ~ 0
  )) %>% 
  filter(remove == 0) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(radiation_start_date =
           coalesce(radiation_start_date,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%

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

  # mutate(treatment_unk = case_when(
  #   str_detect(radiation_start_date, "2300|2301")             ~ "Unknown when and if given 2300",
  #   TRUE                                                      ~ NA_character_
  # )) %>%
   # 3,647 dates created
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




# Bind treatment data----
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






























# Second diagnosis
# Second_dx <- Second_dx %>% 
#   # mutate(deidentified_patient_id = str_to_lower(deidentified_patient_id)) %>% 
#   # filter(!str_detect(deidentified_patient_id, breast_patients))
#   rename(cancer_site = primary_site)
# # write_rds(Second_dx, "Second_dx.rds")
# 
# tbl <- Second_dx %>% 
#   select(cancer_site) %>% 
#   tbl_summary(#sort = list(everything() ~ "frequency")
#   ) %>% as_gt()
# gt::gtsave(tbl, "Cancer site summary.pdf")
# 
# tbl <- Second_dx %>% 
#   select(histology) %>% 
#   tbl_summary(#sort = list(everything() ~ "frequency")
#   ) %>% as_gt()
# gt::gtsave(tbl, "histology summary.pdf")
# 
# 
# Second_dx <- Second_dx %>% 
#   # Create a variable for the first breast diagnose for each patient
#   arrange(deidentified_patient_id, date_of_diagnosis) %>% 
#   mutate(breast_cancer_dx = case_when(
#     str_detect(primary_site, "BREAST")                  ~ date_of_diagnosis
#   )) %>% 
#   select(deidentified_patient_id, breast_cancer_dx, everything()) %>% 
#   arrange(deidentified_patient_id, breast_cancer_dx) %>% 
#   group_by(deidentified_patient_id) %>% 
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
# # distinct(mrn, sample_family_id, sample_id,
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






################################################################################# III ### Merge data
# Get a data with 1 sample per row. Each sample row will have the treatments and demo info
Global_data <- #full_join(Demographic, breast_dna, by = "mrn") %>% 
  # full_join(., breast_info1, by = c("mrn", "date_of_birth")) %>% # I checked the date of birth
  left_join(breast_patients, Treatment, by = "mrn") %>% 
  # full_join(., chemot, by = "mrn") %>% 
  # full_join(., hormonet, by = "mrn") %>% 
  # full_join(., immnunot, by = "mrn") %>% 
  # full_join(., radiot, by = "mrn")
  # Create deidentify IDs
  mutate(rad = "breast_study_") %>%
  group_by(mrn) %>% 
  mutate(id = cur_group_id()) %>%
  ungroup() %>%
  mutate(zero = 6 - nchar(id)) %>%
  mutate(ii = stringi::stri_dup("0", zero)) %>%
  select(c(rad, ii, id, mrn, everything())) %>%
  unite(deidentified_patient_id, rad:id, sep = "") %>% 
  select(c(deidentified_patient_id, mrn, everything(), -zero))

write_rds(Global_data, "Global_data.rds")


# End cleaning
