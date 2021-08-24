# Import Library
library(tidyverse)
library(data.table)
library(VennDiagram)
library(lubridate)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Breast CH")

Demographic <- 
  readxl::read_xlsx(paste0(path, "/raw data/Updated_10R21000238_De-identified_for_Breast_pats_wDNA.xlsx"),
                    sheet = "De-identified_PTE_Demographics") %>% 
  janitor::clean_names()

breast_dna <- 
  readxl::read_xlsx(paste0(path, "/raw data/Updated_10R21000238_De-identified_for_Breast_pats_wDNA.xlsx"),
                    sheet = "UpdatedDeidentified_BioSpecimen") %>% 
  janitor::clean_names()

breast_info <- 
  readxl::read_xlsx(paste0(path, "/raw data/Updated_10R21000238_De-identified_for_Breast_pats_wDNA.xlsx"),
                    sheet = "De-identified_PTE_CR") %>% 
  janitor::clean_names()

Chemot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Updated_10R21000238_De-identified_for_Breast_pats_wDNA.xlsx"),
                    sheet = "De-identified_Treatment_Chemoth") %>% 
  janitor::clean_names()

Hormonet <- 
  readxl::read_xlsx(paste0(path, "/raw data/Updated_10R21000238_De-identified_for_Breast_pats_wDNA.xlsx"),
                    sheet = "De-identified_Treatment_Hormone") %>% 
  janitor::clean_names()

Immnunot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Updated_10R21000238_De-identified_for_Breast_pats_wDNA.xlsx"),
                    sheet = "De-identified_Treatment_Immunot") %>% 
  janitor::clean_names()

Radiot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Updated_10R21000238_De-identified_for_Breast_pats_wDNA.xlsx"),
                    sheet = "De-identified_Treatment_Radiati") %>% 
  janitor::clean_names()



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
  # Make it easier to not have na for future filtering
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
  mutate(chemotherapy_drug = coalesce(chemotherapy_drug, data_unk)) %>% 
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
  # Strat pivot
  select(deidentified_patient_id, chemotherapy_drug, chemotherapy_start_date, chemotherapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, chemotherapy_start_date, chemotherapy_end_date) %>%
  summarise_at(vars(chemotherapy_drug), str_c, collapse = "; ") %>% 
  ungroup()

Chemot <- Chemot %>% 
  select(deidentified_patient_id, treatment_start_date = chemotherapy_start_date, 
         treatment_end_date = chemotherapy_end_date, treatment = chemotherapy_drug) %>% 
  mutate(treatment_type = "chemo")



# pivot wider, use dcast bc better to keep date class
# Chemot <- dcast(setDT(Chemot), deidentified_patient_id ~ rowid(deidentified_patient_id),
#       value.var = c(
#         "chemotherapy_drug",
#         "chemotherapy_start_date",
#         "chemotherapy_end_date"
#       )
# )

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
  mutate(hormone_therapy_drug = coalesce(hormone_therapy_drug, data_unk)) %>% 
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
  # Start pivot
  select(deidentified_patient_id, hormone_therapy_drug, hormone_therapy_start_date, hormone_therapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, hormone_therapy_start_date, hormone_therapy_end_date) %>%
  summarise_at(vars(hormone_therapy_drug), str_c, collapse = "; ") %>% 
  ungroup()

# Hormonet <- dcast(setDT(Hormonet), deidentified_patient_id ~ rowid(deidentified_patient_id),
#                 value.var = c(
#                   "hormone_therapy_drug",
#                   "hormone_therapy_start_date",
#                   "hormone_therapy_end_date"
#                 )
# )

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
  immunotherapy_end_date = as.Date(as.numeric(immunotherapy_end_date), ################################################ Fix bug create a 2300 date, lubridate::origin
                                     origin = "1970-01-01")
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(immunotherapy_start_date = 
           coalesce(immunotherapy_start_date, 
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # Start pivot
  select(deidentified_patient_id, immunotherapy_drug, immunotherapy_start_date, immunotherapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, immunotherapy_start_date, immunotherapy_end_date) %>%
  summarise_at(vars(immunotherapy_drug), str_c, collapse = "; ") %>% 
  ungroup()

# Immnunot <- dcast(setDT(Immnunot), deidentified_patient_id ~ rowid(deidentified_patient_id),
#                   value.var = c(
#                     "immunotherapy_drug",
#                     "immunotherapy_start_date",
#                     "immunotherapy_end_date"
#                   )
# )

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
  mutate(boost_dose_c_gy = coalesce(as.character(boost_dose_c_gy), data_unk)) %>% 
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
                                   origin = "1970-01-01")
  ) %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(radiation_start_date = 
           coalesce(radiation_start_date, 
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>% # 3,647 dates created
  # Start pivot
  select(deidentified_patient_id, boost_dose_c_gy, radiation_start_date, radiation_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, radiation_start_date, radiation_end_date) %>%
  summarise_at(vars(boost_dose_c_gy), str_c, collapse = "; ") %>% 
  ungroup()


# Radiot <- dcast(setDT(Radiot), deidentified_patient_id ~ rowid(deidentified_patient_id),
#                   value.var = c(
#                     "boost_dose_c_gy",
#                     "radiation_start_date",
#                     "radiation_end_date"
#                   )
# )


Radiot <- Radiot %>% 
  select(deidentified_patient_id, treatment_start_date = radiation_start_date, 
         treatment_end_date = radiation_end_date, treatment = boost_dose_c_gy) %>% 
  mutate(treatment_type = "radioT")




# Combine
treatment <- bind_rows(Chemot, Hormonet, Immnunot, Radiot)

# DNA
breast_dna <- breast_dna %>% 
  filter(derived_tissue_type == "Blood",sample_type != "WBC/RBC") %>% 
  mutate(deidentified_patient_id = str_to_lower(deidentified_patient_id)) %>% 
  select(deidentified_patient_id, sample_family_id_sf, sample_id,
         specimen_collection_date) %>% 
  arrange(deidentified_patient_id, specimen_collection_date)


breast_dna <- breast_dna %>% left_join(treatment, by = "deidentified_patient_id") %>% 
  
  mutate(blood_bf_treatment = case_when(
    specimen_collection_date <= treatment_start_date                ~ treatment_type,
    TRUE                                                            ~ "No"
  )) %>% 
  mutate(blood_bf_treatment_30days = case_when(
    specimen_collection_date <= (treatment_start_date + days(30))   ~ treatment_type,
    TRUE                                                            ~ "No"
  )) %>% 
  
  mutate(blood_after_treatment = case_when(
    specimen_collection_date > treatment_start_date                 ~ treatment_type,
    TRUE                                                            ~ "No"
  )) %>% 
  mutate(blood_after_treatment_30days = case_when(
    specimen_collection_date > (treatment_start_date + days(30))   ~ treatment_type,
    TRUE                                                            ~ "No"
  ))




# 
breast_dna <-
  dcast(setDT(breast_dna), deidentified_patient_id+sample_family_id_sf+sample_id+specimen_collection_date ~ rowid(deidentified_patient_id),
        value.var = c("blood_bf_treatment", "blood_bf_treatment_30days",
                      "blood_after_treatment", "blood_after_treatment_30days")) %>% 
  unite(col = c("blood_bf_treatment"), sep = ";", remove = TRUE) %>% 
  unite(col = c("blood_bf_treatment"), sep = ";", remove = TRUE) %>% 
  unite(col = c("blood_bf_treatment"), sep = ";", remove = TRUE) %>% 
  unite(col = c("blood_bf_treatment"), sep = ";", remove = TRUE) 
  
breast_dna <-
  dcast(setDT(breast_dna), deidentified_patient_id ~ rowid(deidentified_patient_id),
        value.var = c("sample_family_id_sf", "sample_id",
                      "specimen_collection_date", ...)) %>%
  
  
  
  
  # purrr::keep(~!all(is.na(.))) %>%
  # select(deidentified_patient_id, ends_with("_1")) %>%
  # `colnames<-`(str_remove(colnames(.), "_1"))
# write_rds(breast_dna, "breast_dna.rds")

breast_info <- 
  breast_info %>% ##################################################### Fix bug
  left_join(., 
            Demographic %>% 
              select(deidentified_patient_id,
                     "date_of_birth"), 
            by = "date_of_birth")








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
