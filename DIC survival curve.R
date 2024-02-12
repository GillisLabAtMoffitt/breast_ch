# Import Library
library(tidyverse)
library(survminer)
library(survival)


################################################################################# I ### Load data
data <- 
  readxl::read_xlsx(paste0(here::here(), "/Data_Supplement Reviewed_07.25.23.xlsx")) %>% 
  mutate(os_event = case_when(
    `Vital status` == "alive"      ~ 0,
    `Vital status` == "dead"      ~ 1
  )) %>% 
  rename(os_time = `Survival time (mo)`)

mysurv <- Surv(time = data$os_time, event = data$os_event)
myplot <- survfit(mysurv~DIC, data = data)

pdf(paste0(here::here(), "/DIC survival plot.pdf"))
ggsurvplot(myplot, data=data,
           title = "Survival Analysis",
           font.main = c(20, "bold", "black"),
           font.x = c(18, "bold", "black"),
           font.y = c(18, "bold", "black"),
           font.legend = c(16, "black"),
           font.tickslab = c(16, "bold", "black"),
           size = 1.5,
           
           xlab = "Time in months",
           legend = "top",
           legend.title = "DIC",
           legend.labs = c("No", "Yes"),
           palette = c("blue", "red"),
           pval = TRUE,
           conf.int = FALSE,
           # Censor
           censor = TRUE,
           # Add risk table
           tables.height = 0.3,
           risk.table.title = "Risk table",
           risk.table = "nrisk_cumevents",
           # cumevents = TRUE,
           risk.table.y.text = FALSE,
           # cumevents.y.text = FALSE,
           risk.table.fontsize = 4,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black"))
) + guides(colour = guide_legend(ncol = 1))
dev.off()
