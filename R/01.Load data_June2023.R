# Import Library
library(tidyverse)
library(data.table)
library(VennDiagram)
library(lubridate)
library(gtsummary)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Breast CH")

Demographic <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "Demographics", na = c("Missing", "Unknown")) %>% 
  janitor::clean_names()

breast_DNA <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_tx_biobabking_report_Yifen_May2023.xlsx"),
                    sheet = "biobanking") %>% 
  janitor::clean_names()

breast_info <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "Cancer Characteristics") %>% 
  janitor::clean_names()

vitals <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "Vital Status ") %>% 
  janitor::clean_names()

treatment <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_tx_biobabking_report_Yifen_May2023.xlsx"),
                    sheet = "clean_all_treatments") %>% 
  janitor::clean_names()

gcsf <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "Medication") %>% 
  janitor::clean_names()

cbc <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "CBC") %>% 
  janitor::clean_names()

breast_marker_1 <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "SSDI_EMR_Labs_for_Breast") %>% 
  janitor::clean_names()

breast_marker_2 <- 
  readxl::read_xlsx(paste0(path, "/raw data/breast_sarcama_report_Yifen_April2023.xlsx"),
                    sheet = "SSF_Site_factor_for_breast") %>% 
  janitor::clean_names()

smoking <- 
  readxl::read_xlsx(paste0(path, "/raw data/smoking_report_for_breast_sarcoma_pts_may72023.xlsx"),
                    sheet = "smoking_cr") %>% 
  janitor::clean_names()


# END Loading

