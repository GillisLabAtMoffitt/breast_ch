library(tidyverse)
library(lubridate)

############################################################ I ### Load data----
path_save <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", 
                      "Breast CH", "sequential_samples_hossein")

sequenced_patient_data <- read_csv(paste0(
  here::here(),
  "/processed data",
  "/Identified breast data with sequenced sequential sample and clinical-CBC-Cytopenia-NADIR_2025-05-16.csv"))

path_raw <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "sequential_samples_hossein")
pfs_chart_rev <- 
  readxl::read_xlsx(
    paste0(path_raw, "/chart_reviewed",
           "/updated PFS date chart review 06-24-25.xlsx")) %>% 
  janitor::clean_names()

her2_chart_rev <- 
  readxl::read_xlsx(
    paste0(path_raw, "/chart_reviewed",
           "/HER2_updated_05.14.2025.xlsx"), na = "NA") %>% 
  janitor::clean_names()

exclude_pfs_other_cancer <- 
  readxl::read_xlsx(
    paste0(path_raw, "/chart_reviewed",
           "/Progression events reviewed_NG_08.14.25.xlsx")) %>% 
  janitor::clean_names()










path_calls <- fs::path("", "Volumes", "Lab_Gillis", "Data", "Breast",
                     "Breast_R01", "CH_calls")

inclusion_data <- readxl::read_xlsx(paste0(path_calls,
                                          "/Breast_final_calls_removeSynonymous_0827_08.12.2025.xlsx"), 
                                   sheet = "Sheet1") %>% 
  janitor::clean_names()






############################################################ II ### Prep new data date----
# PFS
pfs_chart_rev <- pfs_chart_rev %>% 
  `colnames<-`(c("deidentified_patient_id", "mrn",
                 "pfs_time_months_new", "date_of_progression_new", "pfs_date_new")) %>% 
  left_join(., 
            sequenced_patient_data %>% select(mrn, date_of_first_corresponding_treatment),
            by = "mrn") %>% 
  mutate(mrn = as.character(mrn)) %>% 
  mutate(pfs_time_months_new = interval(start = date_of_first_corresponding_treatment, end = pfs_date_new)/
           duration(n = 1, unit = "months")) %>% 
  select(-date_of_first_corresponding_treatment) %>% 
  distinct()
  
  
# HER2
her2_chart_rev <- her2_chart_rev %>% 
  select(deidentified_patient_id, her2) %>% 
  rename(her2_new = her2)
  # left_join(., 
  #           sequenced_patient_data %>% select(mrn, date_of_first_corresponding_treatment),
  #           by = "mrn") %>% 
  # mutate(pfs_time_months_new = interval(start = date_of_first_corresponding_treatment, end = pfs_date_new)/
  #          duration(n = 1, unit = "months")) %>% 
  # select(-date_of_first_corresponding_treatment) %>% 
  # distinct()


############################################################ III ### Replace old date and length/interval variable----
sequenced_patient_data <- sequenced_patient_data %>% 
  mutate(mrn = as.character(mrn)) %>% 
  # update PFS
  left_join(., pfs_chart_rev,
            by = c("mrn", "deidentified_patient_id")) %>% 
  mutate(pfs_time_months = coalesce(pfs_time_months_new, pfs_time_months)) %>% 
  mutate(date_of_progression = coalesce(date_of_progression_new, date_of_progression)) %>% 
  mutate(pfs_date = coalesce(pfs_date_new, pfs_date)) %>% 
  select(-c(pfs_time_months_new, date_of_progression_new, pfs_date_new)) %>% 
  # Update HER2
  left_join(., her2_chart_rev,
            by = "deidentified_patient_id") %>% 
  mutate(her2 = coalesce(her2_new, her2)) %>% 
  select(-her2_new) %>% 
  # Create a new grouped categories variable for Histology
  mutate(
    histology_grouped = case_when(
      #IDC
      histology=="85003 invasive carcinoma of no special type (c50._); 85232 dcis & mixed w/other in situ" ~"Invasive ductal carcinoma",
      histology=="85003 invasive carcinoma of no special type (c50._); 85002 intraductal carcinoma noninfiltrating nos" ~"Invasive ductal carcinoma",
      histology=="85003 invasive carcinoma of no special type (c50._); 85002 intraductal carcinoma noninfiltrating nos; 85003 invasive carcinoma of no special type (c50._)" ~"Invasive ductal carcinoma",
      histology=="85003 invasive carcinoma of no special type (c50._)" ~"Invasive ductal carcinoma",
      histology=="85003 invasive carcinoma of no special type (c50._); 85003 invasive carcinoma of no special type (c50._)" ~"Invasive ductal carcinoma",
      histology=="85003 invasive carcinoma of no special type (c50._); 84013 apocrine adenocarcinoma" ~"Invasive ductal carcinoma",
      histology=="80103 carcinoma nos" ~"Invasive ductal carcinoma",
      histology=="85073 invasive micropapillary carcinoma (c50._)" ~"Invasive ductal carcinoma",
      #ILC
      histology=="85203 lobular carcinoma nos; 85202 lobular carcinoma in situ nos" ~"Invasive lobular carcinoma",
      histology=="85432 paget disease and intraductal carcinoma in situ; 85203 lobular carcinoma nos" ~"Invasive lobular carcinoma",
      histology=="85203 lobular carcinoma nos; 82012 cribriform carcinoma in situ" ~"Invasive lobular carcinoma",
      histology=="85203 lobular carcinoma nos" ~"Invasive lobular carcinoma",
      histology=="85203 lobular carcinoma nos; 85203 lobular carcinoma nos" ~"Invasive lobular carcinoma",
      #Ductal carcinoma in situ (DCIS)
      histology=="82012 cribriform carcinoma in situ" ~"Ductal carcinoma in situ (DCIS)",
      histology=="82302 ductal carcinoma in situ solid type" ~"Ductal carcinoma in situ (DCIS)",
      histology=="85232 dcis & mixed w/other in situ" ~"Ductal carcinoma in situ (DCIS)",
      histology=="85002 intraductal carcinoma noninfiltrating nos" ~"Ductal carcinoma in situ (DCIS)",
      histology=="85002 intraductal carcinoma noninfiltrating nos; 85002 intraductal carcinoma noninfiltrating nos" ~"Ductal carcinoma in situ (DCIS)",
      histology=="85303 inflammatory carcinoma" ~"Ductal carcinoma in situ (DCIS)",
      #Mixed invasive carcinoma (tumors of mixed types)
      histology=="85233 infiltrating duct mixed with other types of carcinoma" ~"Mixed invasive carcinoma (tumors of mixed types)",
      histology=="85223 infiltrating duct and lobular carcinoma" ~"Mixed invasive carcinoma (tumors of mixed types)",
      #Mucinous carcinoma
      histology=="84803 mucinous adenocarcinoma" ~"Special subtypes of invasive carcinoma",
      histology=="84803 mucinous adenocarcinoma; 85003 invasive carcinoma of no special type (c50._)" ~"Special subtypes of invasive carcinoma",
      #Special subtypes of invasive carcinoma
      #Metaplastic carcinoma
      histology=="85753 metaplastic carcinoma nos" ~"Special subtypes of invasive carcinoma",
      #Tubular carcinoma
      histology=="85003 invasive carcinoma of no special type (c50._); 82113 tubular adenocarcinoma" ~"Special subtypes of invasive carcinoma",
      histology=="82113 tubular adenocarcinoma" ~"Special subtypes of invasive carcinoma",
      #Invasive papillary carcinoma
      histology=="80503 papillary carcinoma nos; 85003 invasive carcinoma of no special type (c50._)" ~"Special subtypes of invasive carcinoma"
      
    ), .after = histology
  )


# Exclude PFS that are related to other cancer
sequenced_patient_data <- sequenced_patient_data %>% 
  left_join(., exclude_pfs_other_cancer %>% 
              mutate(mrn = as.character(mrn)),
            by = c("mrn", "progression_type")) %>% 
  
  # recode event 
  mutate(pfs_event = case_when(
    nancy_notes == "Does not count - Different cancer" &
      os_event == 0                                                ~ 0,
    nancy_notes == "Does not count - Different cancer" &
      os_event == 1                                                ~ 1,
    !is.na(progression_type)                                       ~ 1,
    os_event == 1                                                  ~ 1,
    is.na(progression_type)                                        ~ 0
  ), .after = progression_type) %>% 
  unite(progression_type, c(nancy_notes, progression_type), sep = " : ", remove = FALSE, na.rm = TRUE) %>% 
  mutate(progression_type = na_if(progression_type, "")) %>% 

  mutate(new_progression_date = as.Date(new_progression_date)) %>% 
  mutate(date_of_progression = case_when(
    # this will remove the progression date for other cancer
    nancy_notes == "Does not count - Different cancer"           ~ NA_Date_,
    !is.na(date_of_progression)                                  ~ coalesce(new_progression_date, date_of_progression)
  )) %>% 
  mutate(pfs_date = case_when(
    pfs_event == 1                              ~ coalesce(date_of_progression, date_of_last_followup),
    pfs_event == 0                              ~ date_of_last_followup
  )) %>% 
  mutate(pfs_time_months = interval(start = date_of_first_corresponding_treatment, end = pfs_date)/
           duration(n = 1, unit = "months"),
         .after = pfs_event) %>% 

  select(-c(nancy_notes, new_progression_date))


# Add patient inclusion in Mona's analysis

inclusion_data <- inclusion_data %>% 
  select(mrn, patient_included_in_mona_manuscript = included_in_mona_s_manuscript, 
         sample_name_secondary) %>% 
  mutate(mrn = as.character(mrn)) %>% 
  distinct()

sequenced_patient_data <- sequenced_patient_data %>% 
  full_join(., inclusion_data, 
            by = c("mrn", "sample_name_secondary")) %>% 
  select(mrn, deidentified_patient_id, 
         patient_included_in_mona_manuscript, everything()) %>% 
  mutate(samples_excluded = case_when(
    is.na(patient_included_in_mona_manuscript)         ~ "TRUE",
    TRUE                                               ~ "FALSE"
  ), .after = patient_included_in_mona_manuscript) %>% 
  group_by(mrn) %>% 
  fill(patient_included_in_mona_manuscript, .direction = "updown") %>% 
  ungroup()


write_csv(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data_",
                                         today(), ".csv"))
write_csv(sequenced_patient_data, paste0(path_save, "/processed data",
                                         "/Identified breast data_",
                                         today(), ".csv"))
write_rds(sequenced_patient_data, paste0(path_save, "/processed data",
                                         "/Identified breast data_",
                                         today(), ".rds"))
write_rds(sequenced_patient_data, paste0(here::here(), "/processed data",
                                         "/Identified breast data_",
                                         today(), ".rds"))


deids_data <- sequenced_patient_data %>% 
  select(-c(mrn, sample_id, sample_family_id, sample_name_secondary, party_id,
            submitted_ID_for_added_samples, received_ID_for_added_samples,
            Sample_Name_Fastq_file,
            contains("date"),
            blood_bf_chemo_rad : sequenced_samples
  ))

write_csv(deids_data, 
          paste0(path_save, 
                 "/processed data", 
                 "/De-identified breast_", 
                 today(), ".csv"))
write_csv(deids_data, 
          paste0("processed data", 
                 "/De-identified breast_", 
                 today(), ".csv"))

  
  
  
  
  
  
  
  
  