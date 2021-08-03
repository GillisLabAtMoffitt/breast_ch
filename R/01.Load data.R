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
Chemot %>% 
  mutate(chemotherapy_start_date1 = as.Date(chemotherapy_start_date))
  mutate(chemotherapy_start_date1 = as.POSIXct(chemotherapy_start_date), tryFormats = c("%d/%m/%y", "%d/%m/%Y %H:%M:%OS"))











################################################################################# III ### Merge data
Global_data <- full_join(Demographic, breast_dna, by = "deidentified_patient_id")








