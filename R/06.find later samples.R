# Import Library
library(tidyverse)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Breast CH")

sequenced_samples <- 
  read.delim(paste0(path, "/raw data/sample_which_were_sequenced_YH.txt"), 
             header = F)

submitted_samples <- 
  readxl::read_xlsx(paste0(path, "/raw data/samples_submitted_NG.xlsx"),
                    skip = 9) %>% 
  janitor::clean_names() %>% 
  select(sample_name_fastq_file, sample_name_secondary)

blood_patients <- 
  read_rds("/Users/colinccm/Documents/GitHub/Gillis/breast_ch/blood_patients.rds")


################################################################################# II ### Merge data
samples <- sequenced_samples %>% 
  mutate(sample_timing = 
           str_match(
             V1, "(b.*)_(p.*|2.*|3.*)_(chemorad|chemo|radiation|hormone)")[,3]) %>% 
  mutate(sample_treatment_type = 
           str_match(
             V1, "(b.*)_(p.*|2.*|3.*)_(chemorad|chemo|radiation|hormone)")[,4]) %>% 
  mutate(sample_deidentified_id = 
           str_match(
             V1, "(b.*)_(p.*|2.*|3.*)_(chemorad|chemo|radiation|hormone)")[,2]) %>% 
  mutate(id_link_to_sequenced = 
           str_match(
             V1, "(b.*)_(p.*|2.*|3.*)_(chemorad|chemo|radiation|hormone)")[,1]) %>% 
  left_join(., submitted_samples %>% 
              distinct(),
            by= c("id_link_to_sequenced" = "sample_name_fastq_file")) %>% 
  distinct(sample_deidentified_id, id_link_to_sequenced, sample_name_secondary, sample_timing, sample_treatment_type)

rm(sequenced_samples, submitted_samples)

data <- blood_patients %>% 
  mutate(sample_count = sapply(str_split(sample_id, ";"), length)) %>% 
  separate_wider_delim(cols = sample_id, delim = "; ",
                        names = c(paste("sample_id", 1:max(.$sample_count), sep = "_")), 
                        too_few = "align_start", too_many = "drop", 
                        cols_remove = TRUE) %>% 
  pivot_longer(cols = c(starts_with("sample_id_")), 
               names_to = NULL, values_to = "sample_id", 
               values_drop_na = TRUE) %>% 
  select(deidentified_patient_id, mrn, party_id, 
         sample_family_id, sample_id,
         everything(), -sample_count) %>% 
  mutate(sample_id = str_to_upper(sample_id)) %>% 
  full_join(., samples,
            by= c("sample_id" = "sample_name_secondary")) %>% 
  select(mrn, sample_id, id_link_to_sequenced, sample_timing, sample_treatment_type, everything())

data1 <- data %>% 
  mutate(patient_sequenced = case_when(
    !is.na(id_link_to_sequenced)            ~ "Yes"
  ), .before = 4) %>% 
  group_by(mrn) %>% 
  fill(patient_sequenced, .direction = "updown") %>% 
  ungroup() %>% 
  filter(patient_sequenced == "Yes")

# then need to add back the dates















