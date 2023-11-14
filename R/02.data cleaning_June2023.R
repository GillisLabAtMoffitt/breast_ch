################################################################################# II ### Data cleaning

# Breast samples----
table(breast_DNA$sample_type)

breast_dna <- breast_DNA %>% 
  # Select sample type needed for CH detection
  filter(collection_site_anatomic_desc == "Blood", 
         str_detect(sample_type, "Buffy|Genomic|Unprocessed|CD138|MNC$")) %>% 
  select(patient_id, sample_family_id, sample_id,
         specimen_collection_dt) %>%
  # add same sample/same date on the same row
  arrange(patient_id, specimen_collection_dt) %>% 
  # Summarize to have 1 sample/day per row and not 1 row for each aliquot of the same sample 
  group_by(patient_id, sample_family_id, specimen_collection_dt) %>% 
  summarise_at(vars(sample_id), str_c, collapse = "; ") %>%
  ungroup()

write_rds(breast_dna, "breast_dna_06272023.rds")

# ER/PR/HER----
breast_marker_1 <- breast_marker_1 %>% 
  mutate(patient_id = as.character(x1)) %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  select(-c(x1, schema_desc, estrogen_receptor_summary_cd,
            her2_overall_summary_cd, progesterone_receptor_summary_cd)
         ) %>% 
  mutate(ER_PR_status = case_when(
    str_detect(estrogen_receptor_summary_desc, "NEGATIVE") &
      str_detect(progesterone_receptor_summary_desc, "NEGATIVE")        ~ "ER-/PR-",
    str_detect(estrogen_receptor_summary_desc, "NEGATIVE") &
      str_detect(progesterone_receptor_summary_desc, "POSITIVE")        ~ "ER-/PR+",
    str_detect(estrogen_receptor_summary_desc, "POSITIVE") &
      str_detect(progesterone_receptor_summary_desc, "NEGATIVE")        ~ "ER+/PR-",
    str_detect(estrogen_receptor_summary_desc, "POSITIVE") &
      str_detect(progesterone_receptor_summary_desc, "POSITIVE")        ~ "ER+/PR+",
    str_detect(estrogen_receptor_summary_desc, "NEGATIVE")              ~ "ER-",
    str_detect(progesterone_receptor_summary_desc, "NEGATIVE")          ~ "PR-",
    str_detect(estrogen_receptor_summary_desc, "POSITIVE")              ~ "ER+",
    str_detect(progesterone_receptor_summary_desc, "POSITIVE")          ~ "PR+",
  )) %>% 
  mutate(ER_PR_HER_status = case_when(
    ER_PR_status == "ER-/PR-" &
      str_detect(her2_overall_summary_desc, "NEGATIVE")         ~ "ER-/PR-/HER-",
    ER_PR_status == "ER+/PR-" &
      str_detect(her2_overall_summary_desc, "NEGATIVE")         ~ "ER+/PR-/HER-",
    ER_PR_status == "ER-/PR+" &
      str_detect(her2_overall_summary_desc, "NEGATIVE")         ~ "ER-/PR+/HER-",
    ER_PR_status == "ER+/PR+" &
      str_detect(her2_overall_summary_desc, "NEGATIVE")         ~ "ER+/PR+/HER-",
    ER_PR_status == "ER-/PR-" &
      str_detect(her2_overall_summary_desc, "POSITIVE")         ~ "ER-/PR-/HER+",
    ER_PR_status == "ER+/PR-" &
      str_detect(her2_overall_summary_desc, "POSITIVE")         ~ "ER+/PR-/HER+",
    ER_PR_status == "ER-/PR+" &
      str_detect(her2_overall_summary_desc, "POSITIVE")         ~ "ER-/PR+/HER+",
    ER_PR_status == "ER+/PR+" &
      str_detect(her2_overall_summary_desc, "POSITIVE")         ~ "ER+/PR+/HER+",
    ER_PR_status == "ER+" &
      str_detect(her2_overall_summary_desc, "POSITIVE")         ~ "ER+/HER+",
    ER_PR_status == "ER-" &
      str_detect(her2_overall_summary_desc, "POSITIVE")         ~ "ER-/HER+",
    ER_PR_status == "PR+" &
      str_detect(her2_overall_summary_desc, "POSITIVE")         ~ "PR+/HER+",
    ER_PR_status == "PR-" &
      str_detect(her2_overall_summary_desc, "POSITIVE")         ~ "PR-/HER+",
    ER_PR_status == "ER+" &
      str_detect(her2_overall_summary_desc, "NEGATIVE")         ~ "ER+/HER-",
    ER_PR_status == "ER-" &
      str_detect(her2_overall_summary_desc, "NEGATIVE")         ~ "ER-/HER-",
    ER_PR_status == "PR+" &
      str_detect(her2_overall_summary_desc, "NEGATIVE")         ~ "PR+/HER-",
    ER_PR_status == "PR-" &
      str_detect(her2_overall_summary_desc, "NEGATIVE")         ~ "PR-/HER-",
    str_detect(her2_overall_summary_desc, "NEGATIVE")           ~ "HER-",
    str_detect(her2_overall_summary_desc, "POSITIVE")           ~ "HER+"
  )) %>% 
  mutate_at(c("response_to_neoadjuvant_therapy_desc"), ~str_to_sentence(.)) %>% 
  distinct()

breast_marker_2 <- breast_marker_2 %>% 
  rename(ER_results = ssf1_result_desc,
         PR_results = ssf2_result_desc,
         HER_results = ssf9_result_desc) %>% 
  mutate(ER_PR_status = case_when(
    str_detect(ER_results, "Negative/normal") &
      str_detect(PR_results, "Negative/normal")           ~ "ER-/PR-",
    str_detect(ER_results, "Negative/normal") &
      str_detect(PR_results, "Positive/elevated")         ~ "ER-/PR+",
    str_detect(ER_results, "Positive/elevated") &
      str_detect(PR_results, "Negative/normal")           ~ "ER+/PR-",
    str_detect(ER_results, "Positive/elevated") &
      str_detect(PR_results, "Positive/elevated")         ~ "ER+/PR+",
    str_detect(ER_results, "Negative/normal")             ~ "ER-",
    str_detect(PR_results, "Negative/normal")             ~ "PR-",
    str_detect(ER_results, "Positive/elevated")           ~ "ER+",
    str_detect(PR_results, "Positive/elevated")           ~ "PR+",
  )) %>% 
  mutate(ER_PR_HER_status = case_when(
    ER_PR_status == "ER-/PR-" &
      str_detect(HER_results, "Negative/normal")          ~ "ER-/PR-/HER-",
    ER_PR_status == "ER+/PR-" &
      str_detect(HER_results, "Negative/normal")          ~ "ER+/PR-/HER-",
    ER_PR_status == "ER-/PR+" &
      str_detect(HER_results, "Negative/normal")          ~ "ER-/PR+/HER-",
    ER_PR_status == "ER+/PR+" &
      str_detect(HER_results, "Negative/normal")          ~ "ER+/PR+/HER-",
    ER_PR_status == "ER-/PR-" &
      str_detect(HER_results, "Positive/elevated")        ~ "ER-/PR-/HER+",
    ER_PR_status == "ER+/PR-" &
      str_detect(HER_results, "Positive/elevated")        ~ "ER+/PR-/HER+",
    ER_PR_status == "ER-/PR+" &
      str_detect(HER_results, "Positive/elevated")        ~ "ER-/PR+/HER+",
    ER_PR_status == "ER+/PR+" &
      str_detect(HER_results, "Positive/elevated")        ~ "ER+/PR+/HER+",
    ER_PR_status == "ER+" &
      str_detect(HER_results, "Positive/elevated")        ~ "ER+/HER+",
    ER_PR_status == "ER-" &
      str_detect(HER_results, "Positive/elevated")        ~ "ER-/HER+",
    ER_PR_status == "PR+" &
      str_detect(HER_results, "Positive/elevated")        ~ "PR+/HER+",
    ER_PR_status == "PR-" &
      str_detect(HER_results, "Positive/elevated")        ~ "PR-/HER+",
    ER_PR_status == "ER+" &
      str_detect(HER_results, "Negative/normal")          ~ "ER+/HER-",
    ER_PR_status == "ER-" &
      str_detect(HER_results, "Negative/normal")          ~ "ER-/HER-",
    ER_PR_status == "PR+" &
      str_detect(HER_results, "Negative/normal")          ~ "PR+/HER-",
    ER_PR_status == "PR-" &
      str_detect(HER_results, "Negative/normal")          ~ "PR-/HER-",
    str_detect(HER_results, "Negative/normal")            ~ "HER-",
    str_detect(HER_results, "Positive/elevated")          ~ "HER+"
  ))
  
breast_marker <- bind_rows(breast_marker_2, breast_marker_1) %>% 
  distinct(tumor_id, .keep_all = TRUE)
  
rm(breast_marker_1, breast_marker_2)


# Cancer Characteristics----
breast_info <- breast_info %>%
  select(-c(birth_dt : race_cr_src_desc_1, 
            suspense_desc)) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.))) %>% 
  filter(str_detect(primary_site_group_desc, "Breast")) %>%
  distinct(patient_id, dx_dt, .keep_all = TRUE) %>% 
  arrange(patient_id, tumor_id, dx_dt) # I made sure that 1 tumor_id has 1 dx_dt

breast_marker <- breast_marker %>% 
  # Merge with breast info to get date of diagnosis
  left_join(., breast_info %>% 
              select(patient_id, tumor_id, dx_dt),
            by = c("patient_id", "tumor_id"))

# Continue to clean cancer characteristics
breast_info <- breast_info %>% 
  # Summarize to have 1 row per patients
  arrange(patient_id, dx_dt) %>% 
  group_by(patient_id, last_contact_or_death_dt) %>%
  
  summarise_at(vars(dx_dt, tumor_seq_num, 
                    class_of_case_cd, 
                    histology_desc, laterality_desc,
                    stage_clinical_tnm_group_desc, 
                    stage_tnm_cs_mixed_group_desc,
                    metastatic_site_at_diagnosis_desc,
                    type_of_first_recurrence_desc), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  mutate(dx_count = sapply(strsplit(dx_dt, "; "), length), .before = 5) %>% 
  
  separate_wider_delim(dx_dt, 
                       delim = "; ", 
                       names = c(paste("date_of_diagnosis", 1:max(.$dx_count), sep = "")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  separate_wider_delim(stage_tnm_cs_mixed_group_desc, 
                       delim = "; ", 
                       names = c(paste("tnm_stage", 1:max(.$dx_count), sep = "")), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  mutate(across(starts_with("date_of_diagnosis"), ~ as.Date(.))) %>% 
  select(-c(dx_count, last_contact_or_death_dt)) # last_contact_or_death_dt has less rows than in vitals data

breast_info <- breast_info %>% 
  # Merge with markers to have marker status related to the first breast tumor diagnosed
  left_join(., breast_marker,
            by= c("patient_id", "date_of_diagnosis1" = "dx_dt"))

# Demographic----
Demographic <- Demographic %>%
  mutate(mrn = as.character(mrn)) %>% 
  select(patient_id, mrn, birth_dt, gender = cons_gender_derived_desc,
         race = cons_race_derived_desc, 
         ethnicity = cons_ethnicity_derived_desc) %>% 
  mutate(across(where(is.character), ~str_to_sentence(.)))

# Somking
smoking <- smoking %>% 
  mutate(smoking_status = case_when(
    str_detect(derived_tobacco_smoking_status_desc, 
               "Current|Former") ~ "Ever",
    str_detect(derived_tobacco_smoking_status_desc, 
               "Never") ~ "Never",
  )) %>% 
  select(patient_id, smoking_status)


# Vitals----
vitals <- vitals %>% 
  mutate(new_last_date = as.Date(last_contact_or_death_dt, origin = "1899-12-30")) %>% 
  select(patient_id, new_vital_status = cons_derived_vital_status_desc,
         new_last_date)
  

rm(breast_DNA, breast_marker)


################################################################################# III ### Merge data
breast_patients <- breast_dna %>% 
  # Merge with Cancer Char for patients with samples available
  left_join(., breast_info, by = c("patient_id")) %>% 
  # Merge with Demographic for patients with samples available
  left_join(., Demographic, 
            c("patient_id")) %>% 
  left_join(., vitals, 
           c("patient_id")) %>% 
  left_join(., smoking, 
            by= c("patient_id"))
write_rds(breast_patients, "breast_patients_06272023.rds")

breast_patients_id <- paste0(breast_patients$mrn, collapse = "|")


################################################################################# IV ### Clean treatments
# G-CSF----
gcsf <- gcsf %>% 
  distinct(patient_id, drug_catalog_nm) %>% 
  group_by(patient_id) %>%
  summarise_at(vars(drug_catalog_nm), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  mutate(received_gcsf = "Received G-CSF") %>% 
  select(patient_id, received_gcsf, gcsf_type = drug_catalog_nm)


# CBC----
cbc <- cbc %>% 
  filter(lab_nm == "WBC(k/uL)") %>% 
  mutate(neutropenia_at_anytime = case_when(
    lab_result < 1.5                     ~ "Neutropenia (WBC < 1.5 k/uL)"
  )) %>% 
  group_by(patient_id) %>% 
  fill(neutropenia_at_anytime, .direction = "updown") %>% 
  ungroup() %>% 
  distinct(patient_id, .keep_all = TRUE) %>% 
  mutate(neutropenia_at_anytime = case_when(
    !is.na(neutropenia_at_anytime)       ~ neutropenia_at_anytime,
    is.na(neutropenia_at_anytime)        ~ "No"
  )) %>% 
  select(patient_id, neutropenia_at_anytime)


# Chemotherapy----
Chemot <- treatment %>% 
  filter(treatment_type == "CHEMO") %>% 
  # Limit to Breast cancer patients
  # filter(str_detect(patient_id, breast_patients_id)) %>% 
  filter(!is.na(treatment_start_dt) | !is.na(treatment_drug_desc) |
           treatment_rx_summary_desc == "MULTI-AGENT CHEMO")

Chemot1 <- Chemot %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(treatment_start_dt =
           coalesce(treatment_start_dt,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%
  distinct(patient_id, treatment_drug_desc, 
           treatment_start_dt, treatment_end_dt, 
           .keep_all = TRUE) %>% 
  mutate(across(where(is.character), ~str_to_lower(.)))

Chemot2 <- Chemot1 %>%
  arrange(patient_id, treatment_start_dt, treatment_drug_desc) %>% 
  # Combine drugs into regimen
  group_by(patient_id, treatment_start_dt) %>%
  summarise_at(vars(treatment_drug_desc, treatment_end_dt), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(treatment_end_dt, paste("treatment_end_dt", 10:1, sep = ""), 
           sep = "; ", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  mutate(across(starts_with("treatment_end_dt"), ~ as.Date(.))) %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  
  # Remove drugs that are before breast diagnosis
  left_join(., 
            breast_patients %>% 
              distinct(patient_id, .keep_all = TRUE) %>% 
              select(patient_id, date_of_diagnosis1), by= "patient_id") %>% 
  mutate(remove = case_when(
    treatment_start_dt < date_of_diagnosis1             ~ "remove",
    TRUE                                                     ~ NA_character_
  )) %>% 
  filter(is.na(remove))
  
Chemot <- Chemot2 %>% 
  # Combine AC + P as 1 regimen
  group_by(patient_id) %>%
  mutate(linenumber = row_number()) %>%
  
  mutate(AC = case_when(
    treatment_drug_desc == "adriamycin; cyclophosphamide" &
      linenumber == 1                                                 ~ "AC"
  )) %>%
  mutate(AC_start_date = case_when(
    AC == "AC"                                                        ~ treatment_start_dt
  )) %>% 
  mutate(AC_stop_date = case_when(
    AC == "AC"                                                        ~ treatment_end_dt1
  )) %>% 
  mutate(paclitaxel = case_when(
    treatment_drug_desc == "paclitaxel" &
      linenumber == 2                                                 ~ "paclitaxel"
  )) %>%
  mutate(pac_start_date = case_when(
    paclitaxel == "paclitaxel"                                        ~ treatment_start_dt
  )) %>% 
  mutate(pac_end_date = case_when(
    paclitaxel == "paclitaxel"                                        ~ treatment_end_dt1
  )) %>% 
  fill(AC, AC_start_date, AC_stop_date, 
       paclitaxel, pac_start_date, pac_end_date,
       .direction = "updown") %>%
  
  # Combine 5 fu; cyclophosphamide; epirubicin + D as 1 regimen
  mutate(FuCE = case_when(
    treatment_drug_desc == "5 fu; cyclophosphamide; epirubicin" &
      linenumber == 1                                                 ~ "FuCE"
  )) %>%
  mutate(FuCE_start_date = case_when(
    FuCE == "FuCE"                                                    ~ treatment_start_dt
  )) %>% 
  mutate(FuCE_stop_date = case_when(
    FuCE == "FuCE"                                                    ~ treatment_end_dt1
  )) %>% 
  mutate(docetaxel = case_when(
    treatment_drug_desc == "docetaxel" &
      linenumber == 2                                                 ~ "docetaxel"
  )) %>%
  mutate(doce_start_date = case_when(
    docetaxel == "docetaxel"                                          ~ treatment_start_dt
  )) %>% 
  mutate(doce_end_date = case_when(
    docetaxel == "docetaxel"                                          ~ treatment_end_dt1
  )) %>% 
  fill(FuCE, FuCE_start_date, FuCE_stop_date, 
       docetaxel, doce_start_date, doce_end_date,
       .direction = "updown") %>%
  
  #####
  mutate(line_start_gap = as.numeric(treatment_start_dt - lag(treatment_start_dt))) %>% 
  mutate(line_gap = case_when(
    AC == "AC" &
      paclitaxel == "paclitaxel"                                      ~ pac_start_date - AC_stop_date,
    FuCE == "FuCE" &
      docetaxel == "docetaxel"                                        ~ doce_start_date - FuCE_stop_date
  )) %>% 
  fill(line_start_gap, line_gap, .direction = "up") %>% 
  
  ungroup() %>% 
  
  mutate(treatment_drug_desc = case_when(
    (linenumber == 1 |
       linenumber == 2) &
      # line_gap < 90 &
      AC == "AC" &
      paclitaxel == "paclitaxel"                                      ~ "adriamycin; cyclophosphamide; paclitaxel",
    (linenumber == 1 |
       linenumber == 2) &
      FuCE == "FuCE" &
      docetaxel == "docetaxel"                                        ~ "5 fu; cyclophosphamide; epirubicin; docetaxel",
    TRUE                                                              ~ treatment_drug_desc
  )) %>% 
  mutate(treatment_start_dt = case_when(
    (linenumber == 1 |
       linenumber == 2) &
      # line_gap < 90 &
      AC == "AC" &
      paclitaxel == "paclitaxel"                                      ~ AC_start_date,
    (linenumber == 1 |
       linenumber == 2) &
      FuCE == "FuCE" &
      docetaxel == "docetaxel"                                        ~ FuCE_start_date,
    TRUE                                                              ~ treatment_start_dt
  )) %>% 
  mutate(treatment_end_dt1 = case_when(
    (linenumber == 1 |
       linenumber == 2) &
      # line_gap < 90 &
      AC == "AC" &
      paclitaxel == "paclitaxel"                                      ~ pac_end_date,
    (linenumber == 1 |
       linenumber == 2) &
      FuCE == "FuCE" &
      docetaxel == "docetaxel"                                        ~ doce_end_date,
    TRUE                                                              ~ treatment_end_dt1
  )) %>% 
  distinct(patient_id, treatment_drug_desc, treatment_start_dt, treatment_end_dt1, 
           .keep_all = TRUE) %>%  # 705976
  # select(-c(remove, date_of_diagnosis1, AC, AC_start_date, paclitaxel, pac_end_date)) %>% 
  # mutate(linenumber = row_number()) %>% 
  select(patient_id, chemotherapy_drug = treatment_drug_desc, 
         chemotherapy_start_date = treatment_start_dt, 
         chemotherapy_end_date = treatment_end_dt1,
         everything(),
         -c(date_of_diagnosis1, remove))
  

# pivot wider, use dcast bc better to keep date class
chemot <- dcast(setDT(Chemot), patient_id ~ rowid(patient_id),
                value.var = c(
                  "chemotherapy_drug",
                  "chemotherapy_start_date",
                  "chemotherapy_end_date",
                  "line_gap"
                )) %>% 
  select(patient_id, starts_with("chemotherapy_drug"),
         starts_with("chemotherapy_start"),
         starts_with("chemotherapy_end"),
         line_gap_1)

Chemot <- Chemot %>%
  select(patient_id, treatment_start_date = chemotherapy_start_date,
         treatment_end_date = chemotherapy_end_date, treatment = chemotherapy_drug) %>%
  mutate(treatment_type = "chemo")



# Hormonetherapy----
Hormonet <- treatment %>% 
  filter(treatment_type == "HORMONE") %>% 
  # Limit to Breast cancer patients
  # filter(str_detect(patient_id, breast_patients_id)) %>% 
  # Transform number to date and 12:00:00 am character as NA
  filter(!is.na(treatment_start_dt) | !is.na(treatment_drug_desc) |
           treatment_rx_summary_desc == "HORMONE ADMINISTERED" |
           treatment_rx_summary_desc == "RECOMMENDED, UNK IF GIVEN")

Hormonet1 <- Hormonet %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(treatment_start_dt =
           coalesce(treatment_start_dt,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%
  distinct(patient_id, treatment_drug_desc, 
           treatment_start_dt, treatment_end_dt, 
           .keep_all = TRUE) %>% 
  mutate(across(where(is.character), ~str_to_lower(.)))

Hormonet2 <- Hormonet1 %>%
  # Remove drugs not for breast cancer
  mutate(treatment_drug_desc = case_when(
    treatment_drug_desc %in% 
      c("dexamethasone", 
        "synthroid, levothyroxine")                               ~ "non-breast",
    TRUE                                                          ~ treatment_drug_desc
  )) %>% 
  arrange(patient_id, treatment_start_dt, treatment_drug_desc) %>% 
  # Combine drugs into regimen
  group_by(patient_id, treatment_start_dt) %>%
  summarise_at(vars(treatment_drug_desc, treatment_end_dt), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(treatment_end_dt, paste("treatment_end_dt", 10:1, sep = ""), 
           sep = "; ", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  mutate(across(starts_with("treatment_end_dt"), ~ as.Date(.))) %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  
  # Remove drugs that are before breast diagnosis
  left_join(., 
            breast_patients %>% 
              distinct(patient_id, .keep_all = TRUE) %>% 
              select(patient_id, date_of_diagnosis1), by= "patient_id") %>% 
  mutate(remove = case_when(
    treatment_start_dt < date_of_diagnosis1             ~ "remove",
    TRUE                                                     ~ NA_character_
  )) %>% 
  filter(is.na(remove)) %>% 
  select(patient_id, hormone_drug = treatment_drug_desc, 
         hormone_start_date = treatment_start_dt, 
         hormone_end_date = treatment_end_dt1,
         everything(),
         -c(date_of_diagnosis1, remove))

hormonet <- dcast(setDT(Hormonet2), patient_id ~ rowid(patient_id),
                  value.var = c(
                    "hormone_drug",
                    "hormone_start_date",
                    "hormone_end_date"
                  ))

Hormonet <- Hormonet2 %>%
  select(patient_id, treatment_start_date = hormone_start_date,
         treatment_end_date = hormone_end_date, treatment = hormone_drug) %>%
  mutate(treatment_type = "hormone")



# Immnunotherapy----
Immnunot <- treatment %>% 
  filter(treatment_type == "IMMUNO") %>% 
  # Limit to Breast cancer patients
  # filter(str_detect(patient_id, breast_patients_id)) %>% 
  filter(!is.na(treatment_start_dt) | !is.na(treatment_drug_desc) |
           treatment_rx_summary_desc == "RECOMMENDED,UNKN IF GIVEN")

Immnunot1 <- Immnunot %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(treatment_start_dt =
           coalesce(treatment_start_dt,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%
  distinct(patient_id, treatment_drug_desc, 
           treatment_start_dt, treatment_end_dt, 
           .keep_all = TRUE) %>% 
  mutate(across(where(is.character), ~str_to_lower(.)))

Immnunot <- Immnunot1 %>%
  arrange(patient_id, treatment_start_dt, treatment_drug_desc) %>% 
  # Combine drugs into regimen
  group_by(patient_id, treatment_start_dt) %>%
  summarise_at(vars(treatment_drug_desc, treatment_end_dt), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(treatment_end_dt, paste("treatment_end_dt", 10:1, sep = ""), 
           sep = "; ", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  mutate(across(starts_with("treatment_end_dt"), ~ as.Date(.))) %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  
  # Remove drugs that are before breast diagnosis
  left_join(., 
            breast_patients %>% 
              distinct(patient_id, .keep_all = TRUE) %>% 
              select(patient_id, date_of_diagnosis1), by= "patient_id") %>% 
  mutate(remove = case_when(
    treatment_start_dt < date_of_diagnosis1             ~ "remove",
    TRUE                                                     ~ NA_character_
  )) %>% 
  filter(is.na(remove)) %>% 
  select(patient_id, immunotherapy_drug = treatment_drug_desc, 
         immunotherapy_start_date = treatment_start_dt, 
         immunotherapy_end_date = treatment_end_dt1,
         everything(),
         -c(date_of_diagnosis1, remove))


immnunot <- dcast(setDT(Immnunot), patient_id ~ rowid(patient_id),
                  value.var = c(
                    "immunotherapy_drug",
                    "immunotherapy_start_date",
                    "immunotherapy_end_date"
                  ))

Immnunot <- Immnunot %>%
  select(patient_id, treatment_start_date = immunotherapy_start_date,
         treatment_end_date = immunotherapy_end_date, treatment = immunotherapy_drug) %>%
  mutate(treatment_type = "immunoT")



# Radiotion----
Radiot <- treatment %>% 
  filter(treatment_type == "RADIATION") %>% 
  # Limit to Breast cancer patients
  # filter(str_detect(patient_id, breast_patients_id)) %>% 
  filter(!is.na(treatment_start_dt) | !is.na(treatment_end_dt) |
           reason_for_no_treatment_desc == "reason_for_no_treatment_desc")

Radiot1 <- Radiot %>% 
  # Create a really early date to add to NAs date to make sure to exclude these patients
  mutate(treatment_start_dt =
           coalesce(treatment_start_dt,
                    as.Date("1700-01-01", origin = "1899-12-30"))) %>%
  distinct(patient_id, treatment_drug_desc, 
           treatment_start_dt, treatment_end_dt, 
           .keep_all = TRUE) %>% 
  mutate(across(where(is.character), ~str_to_lower(.)))

Radiot2 <- Radiot1 %>% 
  arrange(patient_id, treatment_start_dt, treatment_drug_desc) %>% 
  # Combine drugs into regimen
  group_by(patient_id, treatment_start_dt) %>%
  summarise_at(vars(treatment_drug_desc, treatment_end_dt), 
               str_c, collapse = "; ") %>% 
  ungroup() %>% 
  separate(treatment_end_dt, paste("treatment_end_dt", 10:1, sep = ""), 
           sep = "; ", remove = TRUE, 
           extra = "warn", fill = "left") %>% 
  mutate(across(starts_with("treatment_end_dt"), ~ as.Date(.))) %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  
  # Remove drugs that are before breast diagnosis
  left_join(., 
            breast_patients %>% 
              distinct(patient_id, .keep_all = TRUE) %>% 
              select(patient_id, date_of_diagnosis1), by= "patient_id") %>% 
  mutate(remove = case_when(
    treatment_start_dt < date_of_diagnosis1             ~ "remove",
    TRUE                                                     ~ NA_character_
  )) %>% 
  filter(is.na(remove)) %>% 
  select(patient_id,
         radiation_start_date = treatment_start_dt, 
         radiation_end_date = treatment_end_dt1,
         everything(),
         -c(date_of_diagnosis1, remove))

radiot <- dcast(setDT(Radiot2), patient_id ~ rowid(patient_id),
                value.var = c(
                  "radiation_start_date",
                  "radiation_end_date"
                )
)

Radiot <- Radiot2 %>%
  select(patient_id, treatment_start_date = radiation_start_date,
         treatment_end_date = radiation_end_date) %>%
  mutate(treatment_type = "radioT")




# Bind treatment data----
treatment <- bind_rows(Chemot, Hormonet, Immnunot, Radiot) %>%
  arrange(patient_id, treatment_start_date) %>%
  group_by(patient_id, treatment_type) %>%
  mutate(treatment_line = row_number(patient_id)) %>%
  unite(treatment_line, c(treatment_type, treatment_line), sep = "_", remove = FALSE) %>% 
  ungroup()
write_rds(treatment, "treatment_11142023.rds")

Treatment <- full_join(chemot, hormonet, by = "patient_id") %>% 
  full_join(., immnunot, by = "patient_id") %>% 
  full_join(., radiot, by = "patient_id") #%>% 
  # mutate(had_chemo = case_when(
  #   !is.na(chemotherapy_start_date) ~ "Yes"
  # )) %>% 
  # mutate(had_hormone = case_when(
  #   !is.na(hormone_start_date) ~ "Yes"
  # )) %>% 
  # mutate(had_immuno = case_when(
  #   !is.na(immunotherapy_start_date) ~ "Yes"
  # )) %>% 
  # mutate(had_rad = case_when(
  #   !is.na(radiation_start_date) ~ "Yes"
  # )) %>% 
  # mutate(had_chemo_rad = case_when(
  #   !is.na(chemotherapy_start_date) &
  #     !is.na(radiation_start_date) ~ "Yes"
  # ))

write_rds(Treatment, "Treatment_06272023.rds")






################################################################################# III ### Merge data
# Get a data with 1 sample per row. Each sample row will have the treatments and demo, etc info
Global_data <- 
  left_join(breast_patients, Treatment, by = "patient_id") %>% 
  # left_join(., ERPRHER, by = "mrn") %>% 
  full_join(., gcsf, by = "patient_id") %>% 
  full_join(., cbc, by = "patient_id")


write_rds(Global_data, "Global_data_06272023.rds")

blood_patients <- Global_data %>% 
  # filter to patients who have blood samples
  filter(!is.na(specimen_collection_dt))

write_rds(blood_patients, "blood_patients_06272023.rds")


# End cleaning
