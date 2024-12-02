# setwd("E:/Data/workdata/708279/WOI/iPSYCH2015design/work_temp")
# ==============================================================================
# Title: Data management for study on extended application of case-cohort design
# Author: TW
# Date: Nov 2024
# Description: This script defines samples and weights for primary and
# secondary outcomes until 2015 and 2021 based on datasets including information
# on the iPSYCH design and diagnoses from the National Patient Register
# ==============================================================================

# Content ----------------------------------------------------------------------
# Datamanagement
#   merge studybase and diagnoses (initial datasets generated in Stata)
#   input data 'stam_base2015_all.dta' includes id, iPSYCH case/control status, sex, date om first birth, emigration from DK and death
#   input data 'diag1.dta' includes id, date of first diagnosis of relevant disorders by 2015 or later
#   calculate selection probabilities and weights (Kalbfleish & Lawless, inverse probablitity weights for subcohort, 1 for cases)
#   additional datamanagement
#   sample indicators and weight assignment: 
# the indicator sc_aff identifies cohort sc + aff
# the indicator sc_epi identifies cohort sc + epi (restricted to iPSYCH sample with FU to 2015)
# similarly indicators are generated for other primary and secondary outcomes
# the indicator sc_any identifies all cases in iPSYCH sample (by 2015)

# Load libraries ---------------------------------------------------------------
library(haven)
library(tidyverse)

## Import and merge studybase and diagnoses ------------------------------------
setwd("...")
studybase <- read_dta('stam_base2015_all.dta')
diag <- read_dta('diag1.dta')
studybase <- left_join(studybase, diag, by = "xnr")
rm(diag)

## calculate selection probabilities and weights -------------------------------
# selection probabilities
p1 = 30000/1472762 # numbers taken from iPSYCH papers, Pedersen et al. 2017, Bybjerg-Grauholm et al. 2020. 
p2 = 21000/1657449
sf2015_1981_2005 = 1 - (1 - p1)*(1-p2)
sf2015_1981_2008 = p2 
# inverse probability of sampling weights, by birth cohorts
w81_05 = 1/sf2015_1981_2005
w81_08 = 1/sf2015_1981_2008
# last date of follow-up
studybase$last2015 <- as.Date("2015-12-31")
studybase$last2021 <- as.Date("2021-12-31")
# assign future value to all with 'missing' dates of death to ease handling of missing information
studybase <- studybase %>%
  mutate(death_d = if_else(is.na(death_d),as.Date("3000-12-31"),death_d),
         death_1d = death_d,
         affek2015_1d = affek2015,
         adhd2015_1d = adhd2015)

# Add indicators for secondary outcomes: 
# epilepsy, sud, injuries, TBI, T1DM, asthma, migraine, anxiety, death 
# by 2015 and 2021 nested in the iPSYCH sample, and in the full cohort
disorder <- c("epi","sud10","dm1","inj","tbi","ast","mig","anx", "death")
studybase <- studybase %>% 
  mutate(
    across(all_of(paste0(disorder,"_1d")),
           .fns = list(
             `2015i` = ~if_else((.x <= last2015 & !is.na(.x) & 
                                   (case_any==1 | kontrol2015i==1) & !is.na(kontrol2015i)),1,0),
             `2021i` = ~if_else((.x <= last2021 & !is.na(.x) & 
                                   (case_any==1 | kontrol2015i==1) & !is.na(kontrol2015i)),1,0),
             `2015iall` = ~if_else((.x <= last2015 & !is.na(.x)),1,0),
             `2021iall` = ~if_else((.x <= last2021 & !is.na(.x)),1,0)
           ),
           .names = "{substr(.col,1,nchar(.col)-3)}{.fn}"
    ) 
  )

# Add indicators for primary outcomes up to 2021 
# (for follow-up to 2015, use design variables instead)
disorder <- c("scz10","ssd10","aff10","bip10","asd","adh") 
for (sp in disorder){
  new_var2021 <- paste0(sp,"2021i")
  new_var2021all <- paste0(sp,"2021iall")
  original_var <- paste0(sp,"_1d")
  studybase <- studybase %>%
    mutate(!!new_var2021 :=    ifelse((.data[[original_var]] <= last2021 & 
                                         !is.na(.data[[original_var]]) & 
                                         (case_any==1 | kontrol2015i==1) & 
                                         !is.na(kontrol2015i)),1,0),
           !!new_var2021all := ifelse((.data[[original_var]] <= last2021 & 
                                         !is.na(.data[[original_var]])),1,0))
} 

## sample indicators and weight assignment -------------------------------------
# for full cohort, subcohort and iPSYCH samples
# note that assignment of weights depends on the cases included, 
# i.e. weights != 1 should only be assigned to non-case subcohort members.
# assigning weights 1 to cases and 1/sf to non-case subcohort members.
# weights defined for all iPSYCH subsamples in dataframe. 
# Still need to filter on sample for the different scenarios later.

studybase <- studybase %>%
  mutate(
    full     = 1,
    w_full   = 1,
    sc_any   = ifelse(kontrol2015i==1 | case_any==1, 1, 0),
    w_sc_any = case_when(kontrol2015i==1 & case_any==0 & born1981_2005==1 ~ w81_05,
                         kontrol2015i==1 & case_any==0 & born2006_2008==1 ~ w81_08,
                         case_any==1 ~ 1,
                         TRUE ~ NA_real_),
    sc_orig   = ifelse(kontrol2015i==1 | case_orig==1, 1, 0),
    w_sc_orig = case_when(kontrol2015i==1 & case_orig==0 & born1981_2005==1 ~ w81_05,
                          kontrol2015i==1 & case_orig==0 & born2006_2008==1 ~ w81_08,
                          case_orig==1 ~ 1,
                          TRUE ~ NA_real_),
    sc       = ifelse(kontrol2015i==1, 1, 0),
    w_sc     = case_when(kontrol2015i==1 & born1981_2005==1 ~ w81_05,
                         kontrol2015i==1 & born2006_2008==1 ~ w81_08,
                         TRUE ~ NA_real_)
  )

# Define sets as union of iPSYCH2015 subcohort and primary/secondary outcomes 
# among iPSYCH2015 cases. For primary outcomes use iPSYCH design variables
# and weights corresponding to each scenario
disorder <- c("skizo","skizospek","affek","bipol","autism","adhd","epi","sud10","dm1","inj","tbi","ast","mig","anx","death") # or injuries instead of asthma?
for (sp in disorder){
  new_var <- paste0("sc_",sp)
  new_var_w <- paste0("w_sc_",sp)
  original_var <- paste0(sp,"2015i") # note affek2015i and adhd2015i are iPSYCH2015 design variables
  studybase <- studybase %>%
    mutate(
      !!new_var   := ifelse(kontrol2015i==1 | .data[[original_var]]==1, 1, 0),
      !!new_var_w := case_when(kontrol2015i==1 & !(.data[[original_var]]==1 & case_any==1) & born1981_2005==1 ~ w81_05,
                               kontrol2015i==1 & !(.data[[original_var]]==1 & case_any==1) & born2006_2008==1 ~ w81_08,
                               .data[[original_var]]==1 & case_any==1 ~ 1, # assign weight of 1 only if iPSYCH case
                               TRUE ~ NA_real_) 
    )
}

write.csv(studybase, file = "studybase.csv", row.names = FALSE)