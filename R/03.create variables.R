## USE WHEN GET THE SEQUENCING DATA

# Import library
library(tidyverse)
library(lubridate)

# Load data
blood_patients <- read_rds(paste0(here::here(), "/blood_patients02282023.rds")) %>% 
  mutate(sample_id = str_to_upper(sample_id))

path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Breast CH")
Sequencing <- 
  readxl::read_xlsx(paste0(path, "/raw data/sequencing/samples received EXP1581_add de-ID.xlsx")) %>% 
  janitor::clean_names() %>% 
  mutate(mrn = as.character(mrn))

# Merge to get only the patient with samples
paired_patients <- blood_patients %>% 
  inner_join(Sequencing %>% 
               select(deidentified_patient_id, mrn, 
                      sample_id = sample_lv,
                      sample_type, clinical_event), .,
             by = c("deidentified_patient_id", "mrn",
                    "sample_id"))



# Create var

breast_data_with_blood <- paired_patients %>%
  # create age
  mutate(age_at_diagnosis = round(interval(start = date_of_birth, end = date_of_diagnosis1)/
                                    duration(n = 1, units = "years"), 1)
  ) %>%
  mutate(age_at_sample = round(interval(start = date_of_birth, end = specimen_collection_date)/
                                 duration(n = 1, units = "years"), 1)
  ) %>% 
  mutate(year_at_sample = year(specimen_collection_date)) %>% 
  # mutate(days_from_sample_to_) %>% 
  mutate(gender = "Female") %>% 
  mutate(race_cerner = case_when(
    str_detect(race_cerner, "not")                   ~ NA_character_,
    TRUE                                             ~ race_cerner
  )) %>% 
  mutate(race = coalesce(race_cerner, race_cancer_registry_1, race_derived)) %>% 
  mutate(ethnicity_cerner = case_when(
    str_detect(ethnicity_cerner, "not")              ~ NA_character_,
    TRUE                                             ~ ethnicity_cerner
  )) %>% 
  mutate(ethnicity = coalesce(ethnicity_cerner, csv_ethnicity_cancer_registry, ethnicity_derived))
  
write_rds(breast_data_with_blood, "breast_data_with_blood.rds")


# End create variables
