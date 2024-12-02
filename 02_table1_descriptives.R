# setwd("E:/Data/workdata/708279/WOI/iPSYCH2015design/work_temp")
# ==============================================================================
# Title: Descriptives for table 1
# Author: TW
# Date: Nov 2024
# Description: Crude and upweighted numbers for different case-cohort samples
# ==============================================================================

# Load libraries ---------------------------------------------------------------
library(openxlsx)
library(tidyverse)

# Import data
setwd("...")
studybase <- read_csv("studybase.csv")

# Calculate weights ------------------------------------------------------------
# selection probabilities
p1 = 30000/1472762 # numbers taken from iPSYCH papers, Pedersen et al. 2017, Bybjerg-Grauholm et al. 2020. 
p2 = 21000/1657449
sf2015_1981_2005 = 1 - (1 - p1)*(1-p2)
sf2015_1981_2008 = p2 
# inverse probability of sampling weights, by birth cohorts
w81_05 = 1/sf2015_1981_2005
w81_08 = 1/sf2015_1981_2008
# average weight
w_av = 1/(sf2015_1981_2005*(1472762/1657449) + sf2015_1981_2008*((1657449 - 1472762)/1657449))

# Define samples ---------------------------------------------------------------
sc81_05 <- studybase %>% filter(kontrol2015i==1 & born1981_2005==1) %>% select(xnr)
sc06_08 <- studybase %>% filter(kontrol2015i==1 & born2006_2008==1) %>% select(xnr)
sc <- studybase %>% filter(kontrol2015i==1) %>% select(xnr)

data_sc_any <- studybase %>% filter(sc_any==1) %>% select(xnr, w_sc_any)

# Function to upweight numbers in the pop-representative samples --------------- 
scalar_weighted_n <- function(data,w){
  data_sym <- enquo(data)
  data %>% select(xnr) %>%
    summarise(cohort = !!quo_name(data_sym), 
              n = nrow(data), n_w = n*w
    )
}

# Function to upweight numbers in specific cohort (birth-cohort-specific rows) -
vector_weighted_n <- function(data,subset){
  subset_sym <- enquo(subset)
  subset <- data %>%
    filter(!!subset_sym==1) 
  descr_origset <- subset %>%
    summarise(cohort = !!quo_name(subset_sym), 
              n = nrow(subset), n_w = sum(!!sym(paste0("w_",quo_name(subset_sym))))
    )
  return(descr_origset)
}

# Apply functions to calculate raw and upweighted numbers ----------------------
full_w <- scalar_weighted_n(studybase,1) # note that numbers may slightly differ due to register updates
sc81_05_w <- scalar_weighted_n(sc81_05,w81_05)
sc06_08_w <- scalar_weighted_n(sc06_08,w81_08) # same selection prob as second selection
sc_w_av <- scalar_weighted_n(sc,w_av) # average weights

sc_w <- vector_weighted_n(studybase, sc) # birth cohort specific weights
sc_affek <- vector_weighted_n(studybase, sc_affek) # birth cohort specific weights
sc_epi <- vector_weighted_n(studybase, sc_epi)
sc_any_w <- vector_weighted_n(studybase, sc_any) # birth cohort specific weights
sc_any_w_av <- data_sc_any %>% 
  summarise(cohort="sc_any", n=nrow(data_sc_any), n_w=47657*w_av + 93608*1)

# Combine results
numbers_w <- bind_rows(full_w, sc81_05_w, sc06_08_w, sc_w, sc_w_av, sc_affek, sc_epi, sc_any_w, sc_any_w_av)

# Save results for table 1
write.xlsx(numbers_w,"upweight_tab1.xlsx",rowNames=FALSE, overwrite=TRUE)
