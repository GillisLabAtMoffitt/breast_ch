# Import Library
library(tidyverse)
library(data.table)
library(VennDiagram)
library(lubridate)
library(gtsummary)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Breast CH")

Demographic <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast pts with all blood samples_01.27.22.xlsx"),
                    sheet = "PTE Demographics") %>% 
  janitor::clean_names()

breast_DNA <- # Updated_10R21000238_De-identified_for_Breast_pats_wDNA.xlsx
  readxl::read_xlsx(paste0(path, "/raw data/Breast pts with all blood samples_01.27.22.xlsx"),
                    sheet = "All Blood Specimens") %>% 
  janitor::clean_names()

breast_info <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast pts with all blood samples_01.27.22.xlsx"),
                    sheet = "PTE Cancer Characteristics") %>% 
  janitor::clean_names()

Chemot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast pts with all blood samples_01.27.22.xlsx"),
                    sheet = "Treatment_Chemotherapy") %>% 
  janitor::clean_names()

Hormonet <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast pts with all blood samples_01.27.22.xlsx"),
                    sheet = "Treatment_Hormone") %>% 
  janitor::clean_names()

Immnunot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast pts with all blood samples_01.27.22.xlsx"),
                    sheet = "Treatment_Immuno") %>% 
  janitor::clean_names()

Radiot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast pts with all blood samples_01.27.22.xlsx"),
                    sheet = "Treatment_Radiation") %>% 
  janitor::clean_names()

ERPRHER <- 
  readxl::read_xlsx(paste0(path, "/raw data/Breast pts with all blood samples_01.27.22.xlsx"),
                    sheet = "SS_HER2") %>% 
  janitor::clean_names()

# Second_dx <- read_csv(paste0(path, "/raw data/Data_export_for_De-identified_PTE_CR_14-Sep-2021.csv")) %>% 
#   janitor::clean_names() %>% 
#   select(-group_name)

