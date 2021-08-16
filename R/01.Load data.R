# Import Library
library(tidyverse)
library(data.table)
library(VennDiagram)


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
# Domo
Demographic <- Demographic %>% 
  mutate(across(where(is.character), ~str_to_lower(.)))
  

# DNA
breast_dna <- breast_dna %>% 
  filter(derived_tissue_type == "Blood",sample_type != "WBC/RBC") %>% 
  mutate(deidentified_patient_id = str_to_lower(deidentified_patient_id)) %>% 
  select(deidentified_patient_id, sample_family_id_sf, sample_id,
         specimen_collection_date) %>% 
  arrange(deidentified_patient_id, specimen_collection_date)
# Select the earliest sample
breast_dna <-
  dcast(setDT(breast_dna), deidentified_patient_id ~ rowid(deidentified_patient_id), 
        value.var = c("sample_family_id_sf", "sample_id", 
                      "specimen_collection_date")) %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  select(deidentified_patient_id, ends_with("_1")) %>% 
  `colnames<-`(str_remove(colnames(.), "_1"))
  
# write_rds(breast_dna, "breast_dna.rds")

# Chemot
Chemot <- Chemot %>% 
  filter(chemotherapy_type == "CHEMO NOS" |
           chemotherapy_type == "CONTRAINDICATED" |
           chemotherapy_type == "MULTI-AGENT CHEMO" |
           chemotherapy_type == "NONE, NOT PLANNED" |
           chemotherapy_type == "SINGLE-AGENT CHEMO") %>% 
  
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  mutate(chemotherapy_start_date = case_when(
    str_detect(chemotherapy_start_date, "12:00:00 AM")       ~ NA_character_,
    TRUE                                                     ~ chemotherapy_start_date
  ), 
  chemotherapy_start_date = as.Date(as.numeric(chemotherapy_start_date), 
                                                origin = "1899-12-30")
  ) %>% 
  mutate(chemotherapy_end_date = case_when(
    str_detect(chemotherapy_end_date, "12:00:00 AM")         ~ NA_character_,
    TRUE                                                     ~ chemotherapy_end_date
  ), 
  chemotherapy_end_date = as.Date(as.numeric(chemotherapy_end_date), 
                                     origin = "1899-12-30")
  ) %>% 
  select(deidentified_patient_id, chemotherapy_drug, chemotherapy_start_date, chemotherapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, chemotherapy_start_date, chemotherapy_end_date) %>%
  summarise_at(vars(chemotherapy_drug), str_c, collapse = "; ") %>% 
  group_by(deidentified_patient_id) %>%
  summarise_at(vars(chemotherapy_drug, chemotherapy_start_date, chemotherapy_end_date), str_c, collapse = "; ") %>% 
  separate(col = chemotherapy_drug, paste("chemotherapy_drug", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  separate(col = chemotherapy_start_date, paste("chemotherapy_start_date_", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  separate(col = chemotherapy_end_date, paste("chemotherapy_end_date_", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  purrr::keep(~!all(is.na(.)))

# library(lubridate)

# Hormonet
Hormonet <- Hormonet %>% 
  filter(hormone_therapy_type == "HORMONE ADMINISTERED") %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  mutate(hormone_therapy_start_date = case_when(
    str_detect(hormone_therapy_start_date, "12:00:00 AM")    ~ NA_character_,
    TRUE                                                     ~ hormone_therapy_start_date
  ), 
  hormone_therapy_start_date = as.Date(as.numeric(hormone_therapy_start_date), 
                                    origin = "1899-12-30")
  ) %>% 
  mutate(hormone_therapy_end_date = case_when(
    str_detect(hormone_therapy_end_date, "12:00:00 AM")      ~ NA_character_,
    TRUE                                                     ~ hormone_therapy_end_date
  ), 
  hormone_therapy_end_date = as.Date(as.numeric(hormone_therapy_end_date), 
                                  origin = "1899-12-30")
  ) %>% 
  select(deidentified_patient_id, hormone_therapy_drug, hormone_therapy_start_date, hormone_therapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, hormone_therapy_start_date, hormone_therapy_end_date) %>%
  summarise_at(vars(hormone_therapy_drug), str_c, collapse = "; ") %>% 
  group_by(deidentified_patient_id) %>%
  summarise_at(vars(hormone_therapy_drug, hormone_therapy_start_date, hormone_therapy_end_date), str_c, collapse = "; ") %>% 
  separate(col = hormone_therapy_drug, paste("hormone_therapy_drug", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  separate(col = hormone_therapy_start_date, paste("hormone_therapy_start_date_", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  separate(col = hormone_therapy_end_date, paste("hormone_therapy_end_date_", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  purrr::keep(~!all(is.na(.)))

# Immnunot
Immnunot <- Immnunot %>% 
  filter(immunotherapy_type == "IMMUNO ADMINISTERED") %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  mutate(immunotherapy_start_date = case_when(
    str_detect(immunotherapy_start_date, "12:00:00 AM")       ~ NA_character_,
    TRUE                                                      ~ immunotherapy_start_date
  ), 
  immunotherapy_start_date = as.Date(as.numeric(immunotherapy_start_date), 
                                       origin = "1899-12-30")
  ) %>% 
  mutate(immunotherapy_end_date = case_when(
    str_detect(immunotherapy_end_date, "12:00:00 AM")       ~ NA_character_,
    TRUE                                                      ~ immunotherapy_end_date
  ), 
  immunotherapy_end_date = as.Date(as.numeric(immunotherapy_end_date), 
                                     origin = "1899-12-30")
  ) %>% 
  select(deidentified_patient_id, immunotherapy_start_date, immunotherapy_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id) %>%
  summarise_at(vars(immunotherapy_start_date, immunotherapy_end_date), str_c, collapse = "; ") %>% 
  separate(col = immunotherapy_start_date, paste("immunotherapy_start_date_", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  separate(col = immunotherapy_end_date, paste("immunotherapy_end_date_", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  purrr::keep(~!all(is.na(.)))

# Radiot
Radiot <- Radiot %>% 
  filter(reason_for_no_radiation == "RAD THERAPY PERFORMED"| 
         is.na(reason_for_no_radiation)) %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  mutate(immunotherapy_start_date = case_when(
    str_detect(radiation_start_date, "12:00:00 AM")       ~ NA_character_,
    TRUE                                                      ~ radiation_start_date
  ), 
  radiation_start_date = as.Date(as.numeric(radiation_start_date), 
                                     origin = "1899-12-30")
  ) %>% 
  mutate(radiation_end_date = case_when(
    str_detect(radiation_end_date, "12:00:00 AM")       ~ NA_character_,
    TRUE                                                      ~ radiation_end_date
  ), 
  radiation_end_date = as.Date(as.numeric(radiation_end_date), 
                                   origin = "1899-12-30")
  ) %>% 
  select(deidentified_patient_id, radiation_start_date, radiation_end_date) %>% 
  distinct() %>% 
  group_by(deidentified_patient_id, radiation_start_date) %>%
  summarise_at(vars(radiation_end_date), str_c, collapse = "; ") %>% 
  separate(col = radiation_end_date, paste("radiation_end_date_", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  group_by(deidentified_patient_id) %>%
  summarise_at(vars(radiation_start_date, radiation_end_date_1), str_c, collapse = "; ") %>% 
  separate(col = radiation_start_date, paste("radiation_start_date_", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  separate(col = radiation_end_date_1, paste("regimen_rad_end_date_", 1:5, sep=""), sep = "; |;", extra = "warn", 
           fill = "right") %>% 
  purrr::keep(~!all(is.na(.)))

################################################################################# III ### Merge data
Global_data <- full_join(Demographic, breast_dna, by = "deidentified_patient_id") %>% 
  full_join(., Chemot, by = "deidentified_patient_id") %>% 
  full_join(., Hormonet, by = "deidentified_patient_id") %>% 
  full_join(., Immnunot, by = "deidentified_patient_id") %>% 
  full_join(., Radiot, by = "deidentified_patient_id")

write_rds(Global_data, "Global_data.rds")




