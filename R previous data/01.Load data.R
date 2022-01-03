# Import Library
library(tidyverse)
library(data.table)
library(VennDiagram)
library(lubridate)
library(gtsummary)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Breast CH")

Demographic <- 
  readxl::read_xlsx(paste0(path, "/raw data/10R21000238_De-identified_for_Breast_pats_wDNAsample.xlsx"),
                    sheet = "De-identified_PTE_Demographics") %>% 
  janitor::clean_names()

breast_DNA <- 
  readxl::read_xlsx(paste0(path, "/raw data/Updated_10R21000238_De-identified_for_Breast_pats_wDNA.xlsx"),
                    sheet = "UpdatedDeidentified_BioSpecimen") %>% 
  janitor::clean_names()

breast_info <- 
  readxl::read_xlsx(paste0(path, "/raw data/10R21000238_De-identified_for_Breast_pats_wDNAsample.xlsx"),
                    sheet = "De-identified_PTE_CR_V2") %>% 
  janitor::clean_names()

Chemot <- 
  readxl::read_xlsx(paste0(path, "/raw data/10R21000238_De-identified_for_Breast_pats_wDNAsample.xlsx"),
                    sheet = "De-identified_Treatment_Chemoth") %>% 
  janitor::clean_names()

Hormonet <- 
  readxl::read_xlsx(paste0(path, "/raw data/10R21000238_De-identified_for_Breast_pats_wDNAsample.xlsx"),
                    sheet = "De-identified_Treatment_Hormone") %>% 
  janitor::clean_names()

Immnunot <- 
  readxl::read_xlsx(paste0(path, "/raw data/10R21000238_De-identified_for_Breast_pats_wDNAsample.xlsx"),
                    sheet = "De-identified_Treatment_Immunot") %>% 
  janitor::clean_names()

Radiot <- 
  readxl::read_xlsx(paste0(path, "/raw data/10R21000238_De-identified_for_Breast_pats_wDNAsample.xlsx"),
                    sheet = "De-identified_Treatment_Radiati") %>% 
  janitor::clean_names()

Second_dx <- read_csv(paste0(path, "/raw data/Data_export_for_De-identified_PTE_CR_14-Sep-2021.csv")) %>% 
  janitor::clean_names() %>% 
  select(-group_name)

