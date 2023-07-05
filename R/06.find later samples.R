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
  # extract all info from sequenced sample
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

# write_rds(samples, "samples sequenced as of june 2023.rds")
rm(sequenced_samples, submitted_samples)

data <- blood_patients %>% 
  mutate(specimen_collection_date = as.Date(specimen_collection_date)) %>% 
  # Separate back sample ids to pivot longer
  # Need 1 id per row to merge with sample sequenced
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
  # full join with sample sequenced to keep the later samples
  full_join(., samples,
            by= c("sample_id" = "sample_name_secondary")) %>% 
  select(mrn, sample_id, id_link_to_sequenced, 
         sample_timing, sample_treatment_type, 
         specimen_collection_date, everything())

data1 <- data %>% 
  arrange(mrn, specimen_collection_date, id_link_to_sequenced) %>% 
  # Keep only the patients who were sequenced
  mutate(patient_sequenced = case_when(
    !is.na(id_link_to_sequenced)            ~ "Yes"
  ), .before = 4) %>% 
  group_by(mrn) %>% 
  fill(patient_sequenced, .direction = "updown") %>% 
  ungroup() %>% 
  filter(patient_sequenced == "Yes") %>% 
  # Create a temp variable to help keep only the sample which were sequenced
  # and not the one at the same date but not sequenced
  # but also keep all of the not sequenced
  mutate(temp_rank_number = case_when(
    is.na(id_link_to_sequenced)             ~ dense_rank(mrn)
  ), .before = 3) %>% 
  group_by(mrn, specimen_collection_date) %>% 
  fill(temp_rank_number, specimen_collection_date, .direction = "updown") %>% 
  ungroup() %>%
  # Now filter
  distinct(mrn, temp_rank_number, specimen_collection_date, .keep_all = TRUE) %>% 
  arrange(mrn, specimen_collection_date) %>% 
  select(-temp_rank_number)

data2 <- data1 %>% 
  mutate(is_extra_sample = case_when(
    is.na(id_link_to_sequenced)             ~ "Yes",
    !is.na(id_link_to_sequenced)            ~ "No"
  ), .before = 3) %>% 
  # keep only the later sample
  group_by(mrn, is_extra_sample) %>% 
  mutate(last_sample_seq_date = case_when(
    !is.na(id_link_to_sequenced)             ~ last(specimen_collection_date)
  ), .before = 3) %>% 
  group_by(mrn) %>% 
  fill(last_sample_seq_date, .direction = "updown") %>% 
  ungroup() %>% 
  
  
  
  mutate(has_a_later_sample = case_when(
    is.na(id_link_to_sequenced) &
      specimen_collection_date > last_sample_seq_date    ~ "Yes"
  ), .before = 3) %>% 
  group_by(mrn) %>% 
  fill(has_a_later_sample, .direction = "updown") %>% 
  ungroup()




data3 <- data2 %>% 
  # mutate(has_notsequenced_sample = case_when(
  #   is.na(id_link_to_sequenced)             ~ "Yes"
  # )) %>% 
  # group_by(mrn) %>% 
  # fill(has_notsequenced_sample, sample_treatment_type,
  #      .direction = "updown") %>% 
  # ungroup() %>% 
  # filter(has_notsequenced_sample == "Yes") %>% 
  # select(-temp_rank_number)
  filter(has_a_later_sample == "Yes") %>% 
  select(-c(last_sample_seq_date, has_a_later_sample, patient_sequenced))

nrow(data3 %>% distinct(mrn))

data4 <- data3 %>% 
  arrange(mrn, specimen_collection_date) %>% 
  group_by(mrn) %>% 
  fill(sample_treatment_type, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(time_samples_from_treatment = case_when(
    sample_treatment_type == "chemo"       ~ 
      (interval(start = chemotherapy_start_date_1, end = specimen_collection_date) /
         duration(n = 1, units = "days")),
    sample_treatment_type == "radiation"       ~ 
      (interval(start = radiation_start_date_1, end = specimen_collection_date) /
         duration(n = 1, units = "days")),
    sample_treatment_type == "chemorad"       ~ 
      (interval(start = radiation_start_date_1, end = specimen_collection_date) /
         duration(n = 1, units = "days")),
    sample_treatment_type == "hormone"       ~ 
      (interval(start = hormone_therapy_start_date_1, end = specimen_collection_date) /
         duration(n = 1, units = "days"))
  ), .before = 3) %>% 
  filter(time_samples_from_treatment > 0)

nrow(data4 %>% distinct(mrn))

write_csv(data4, "list of patients with later samples.csv")










