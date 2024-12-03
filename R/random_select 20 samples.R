# Import Library
library(tidyverse)
library(here)
library(lubridate)


################################################################################# I ### Load data
blood_patients <- 
  read_rds(paste0(here(), "/blood_patients.rds"))

available1 <- 
  readxl::read_xlsx(paste0(here(), "/EXP1581DNA - Dr. Gillis 21545 Datasheet_w locations.xlsx"),
                    skip = 9, n_max = 74
                    ) %>% 
  janitor::clean_names()

available2 <- 
  readxl::read_xlsx(paste0(here(), "/EXP1581DNA - Dr. Gillis 21545 Datasheet_w locations.xlsx"),
                    sheet = "DNA Extraction Part B-WB",
                    skip = 9#, n_max = 74
  ) %>% 
  janitor::clean_names() %>% 
  rename(aliquot_final_qty_ug = aliquot_final_quantity_ug,
         aliquot_volume_needed_from_master_u_l = volume_needed_for_1ug_u_l) %>% 
  mutate(aliquot_volume_needed_from_master_u_l = as.numeric(aliquot_volume_needed_from_master_u_l))


################################################################################# II ### Merge
available <- bind_rows(available1, available2, .id = "excel_tab") %>% 
  mutate(excel_tab = case_when(
    excel_tab == 1    ~ "first tab",
    excel_tab == 2    ~ "second tab",
  )) %>% 
  select(-c("sample_run_number", "lv_alias_dna_id", "surgical_no",
            "study")) %>% 
  mutate(sample_lv = str_to_lower(sample_lv))

blood_patients <- blood_patients %>% 
  # group_by(mrn, party_id, sample_family_id, specimen_collection_date) %>% 
  separate(col = sample_id, into = paste("sample_id", 1:30, sep="_"), sep = "; ") %>% 
  purrr::keep(~!all(is.na(.))) %>%
# separate(col = drug_name_, paste("drug_name_", 1:7, sep=""), sep = ";", extra = "warn", # Just for 1 row
#          fill = "right")
  arrange(mrn, specimen_collection_date) %>% 
  pivot_longer(cols = starts_with("sample_id"),
               names_to = "sample_sequence", 
               values_to = "sample_id", 
               values_drop_na = TRUE) %>% 
  select("mrn", "deidentified_patient_id", "party_id",
         "sample_family_id", "sample_sequence", "sample_id",
         "specimen_collection_date", everything())

breast_data <- left_join(available #%>% 
                           # select(-c("sample_run_number", "lv_alias_dna_id", "surgical_no",
                           #           "study")) %>% 
                           # mutate(sample_lv = str_to_lower(sample_lv))
                         , 
                         blood_patients%>% 
                           mutate(sample_family_id = str_to_upper(sample_family_id)),
                         by = c("mrn", "sample_family_id", "sample_lv" = "sample_id"))


breast_data <- breast_data %>% 
  group_by(mrn) %>% 
  mutate(number_samples = n()) %>%
  arrange(desc(number_samples), mrn, specimen_collection_date) %>% 
  # distinct(mrn, clinical_event, .keep_all = TRUE) %>% 
  mutate(specimen_sequence = row_number()) %>%
  group_by(mrn) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  group_by(mrn, clinical_event) %>% 
  mutate(specimen_sequence_1patient_have_2_sample_same_day = row_number()) %>%
  ungroup() %>%
  mutate(age_at_sample = round(interval(start = date_of_birth, end = specimen_collection_date)/
                                 duration(n = 1, units = "years"), 1)
  ) %>% 
  select("mrn", "deidentified_patient_id", "party_id",
         age_at_diagnosis, age_at_sample,
         clinical_event, chemotherapy_drug_1,
         specimen_sequence_1patient_have_2_sample_same_day,
         sample_id = "sample_lv", "sample_family_id",
         "specimen_collection_date",
         "immediate_parent_sample_id",
         "dna_master_lv_id",
         "dna_aliquot_lv_id",
         "original_parent_sample_id", everything(),
         -c(number_samples, specimen_sequence, sample_sequence, n)
         )

write_csv(breast_data, paste0(here(), "/sample lists/list of samples available in the bank.csv"))

breast_data <- breast_data %>% 
  filter(!str_detect(clinical_event, "hormone"))

write_csv(breast_data, paste0(here(), "/sample lists/list of samples available in the bank hormone removed.csv"))


