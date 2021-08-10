# Import Library
library(tidyverse)
library(data.table)

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
# DNA
breast_dna <- breast_dna %>% 
  janitor::clean_names() %>% 
  filter(derived_tissue_type == "Blood",sample_type != "WBC/RBC") %>% 
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
  
write_rds(breast_dna, "breast_dna.rds")

# Chemot
a <- Chemot %>% 
  mutate(across(where(is.character), ~str_to_lower(.))) %>% 
  mutate(chemotherapy_start_date = case_when(
    str_detect(chemotherapy_start_date, "12:00:00 AM")    ~ NA_character_,
    TRUE                                                  ~ chemotherapy_start_date
  ), 
  chemotherapy_start_date = as.Date(as.numeric(chemotherapy_start_date), 
                                                origin = "1899-12-30")
  ) %>% 
  mutate(chemotherapy_end_date = case_when(
    str_detect(chemotherapy_end_date, "12:00:00 AM")    ~ NA_character_,
    TRUE                                                  ~ chemotherapy_end_date
  ), 
  chemotherapy_end_date = as.Date(as.numeric(chemotherapy_end_date), 
                                     origin = "1899-12-30")
  ) %>% 
  group_by(deidentified_patient_id, chemotherapy_start_date) %>%
  summarise_at(vars(chemotherapy_drug, chemotherapy_end_date), str_c, collapse = "; ")
  # summarise_at(vars(chemotherapy_drug, chemotherapy_end_date), c(paste), collapse = "; ") %>% 
  # separate(chemotherapy_end_date, "chemotherapy_end_date", sep = "; ", extra = "drop")
  


library(lubridate)



################################################################################# III ### Merge data
Global_data <- full_join(Demographic, breast_dna, by = "deidentified_patient_id")








