# sanity check

# Look at the end date before start date

wrong_date <- blood_patients2 %>% 
  filter(!if_all(starts_with("seqsample_id"), is.na)) %>%
  # select(starts_with("seqsample_id"))
  mutate(wrong_chemo_date = case_when(
    chemotherapy_start_date_1 > chemotherapy_end_date1_1      ~ "wrong"
  )) %>% 
  mutate(wrong_hormone_date = case_when(
    hormone_therapy_start_date_1 > hormone_therapy_end_date1_1      ~ "wrong"
  ))

write_csv(wrong_date %>% 
            filter(wrong_chemo_date == "wrong") %>% 
            select(mrn, chemotherapy_start_date_1, chemotherapy_end_date1_1, chemotherapy_drug_1), 
          "Wrong chemotherapy dates.csv")

write_csv(wrong_date %>% 
            filter(wrong_hormone_date == "wrong") %>% 
            select(mrn, hormone_therapy_start_date_1, hormone_therapy_end_date1_1, hormone_therapy_drug_1), 
          "Wrong hormonetherapy dates.csv")

