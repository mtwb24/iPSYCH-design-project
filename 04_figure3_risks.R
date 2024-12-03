# ==============================================================================
# Title: Smoothed risks across age (Figure 3)
# Author: TW
# Date: Nov 2024
# Description: Cumulative incidence estimates (1-KM) across age for affective
# disorder (primary outcome) and epilepsy (secondary outcome) 
# using GAM (generalized additive models) for smoothing
# Across four samples and with follow-up by calendar years 2015 and 2021
# Estimates calculated for ages 18, 25, 30, (40), and plotted over age.
# ==============================================================================

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(broom)
library(survival)
library(Epi)
library(mgcv) 
library(openxlsx)
library(ggplot2)
library(lemon)

# Import data
studybase <- read_csv("studybase.csv")

## time to event - function ----------------------------------------------------
# select variables, define survival data, indicator and date of event/censoring
studybase_tte <- function(data,eventdvar, age_tstart, lastfu, emidate, sample){
  data %>%
    select(pnr,fdato.x,{{eventdvar}}, death_d,{{emidate}}, born1981_2005, born2006_2008, {{lastfu}},kontrol2015i, affek2015i, 
           case_any, full, w_full, sc_any, w_sc_any, paste0(sample), paste0("w_",sample), sc, w_sc, kqn.x) %>%  
    mutate(tstop = pmin({{eventdvar}}, {{lastfu}}, death_d, {{emidate}}, 
                        na.rm = TRUE)) %>%
    mutate(event = ifelse({{eventdvar}}>tstop | is.na({{eventdvar}}),0,1),
           age_tstart = age_tstart,
           age_tstop = (unclass(tstop) - unclass(fdato.x) ) /365.24,
           id = pnr,
           death = death_d == tstop
           ) 
}

# different data frames for different events of interest (and follow-up)
# affective disorder
studybase_affek2015 <- studybase_tte(studybase,aff_1d,age_tstart=1,last2015,emi10y_d,sample="sc_affek")
studybase_affek2021 <- studybase_tte(studybase,aff_1d, age_tstart=1,last2021,emi10y_d,sample="sc_affek")
# epilepsy
studybase_epi2015 <- studybase_tte(studybase,epi_1d,age_tstart=1,last2015,emi10y_d,sample="sc_epi")
studybase_epi2021 <- studybase_tte(studybase,epi_1d, age_tstart=1,last2021,emi10y_d,sample="sc_epi")

# Fit 1-KM and store estimates in dataframes
years <- c(2015, 2021)
outcomes <- c("affek", "epi")
models <- list()
for (year in years) {
  for (out in outcomes) {
    survfit_full <- survfit(Surv(age_tstop, event) ~ 1, data = get(paste0("studybase_", out, year)))
    assign(paste0("survfit_full_", out, year), survfit_full)
    assign(paste0("survest_full_", out, year),survfit_full %>% 
             broom::tidy() %>% 
             select(time, estimate, conf.high, conf.low, std.error) %>%  
             rename(p = estimate))
    survfit_sc <- survfit(Surv(age_tstop, event) ~ 1, data = get(paste0("studybase_", out, year)), 
                          subset = sc == 1)
    assign(paste0("survfit_sc_", out, year), survfit_sc)
    assign(paste0("survest_sc_", out, year), survfit_sc %>% 
             broom::tidy() %>% 
             select(time, estimate, conf.high, conf.low, std.error) %>%  
             rename(p = estimate))
    survfit_sc_dis <- survfit(Surv(age_tstop, event) ~ 1, data = get(paste0("studybase_", out, year)), 
                          subset = get(paste0("sc_", out)) == 1, weights = get(paste0("w_sc_", out)), cluster = pnr)
    assign(paste0("survfit_sc_", out, "_", out, year), survfit_sc_dis)
    assign(paste0("survest_sc_", out, "_", out, year), survfit_sc_dis %>% 
             broom::tidy() %>% 
             select(time, estimate, conf.high, conf.low, std.error) %>%  
             rename(p = estimate))
    survfit_sc_any <- survfit(Surv(age_tstop, event) ~ 1, data = get(paste0("studybase_", out, year)), 
                          subset = sc_any == 1, weights = w_sc_any, cluster = pnr)
    assign(paste0("survfit_sc_any_", out, year), survfit_sc_any)
    assign(paste0("survest_sc_any_", out, year), survfit_sc_any %>% 
             broom::tidy() %>% 
             select(time, estimate, conf.high, conf.low, std.error) %>%  
             rename(p = estimate))
  }
}

# Risk estimates at specific time points (ages)

# function for storing estimates in data frame
res_times <- function(survfitdata,out, year,time, scenario){
  data.frame(time=time, 
             out=out,
             year=year,
             cuminc = round((1 - summary(survfitdata, time=time)$surv)*100,2), 
             cuminc_low = round((1 - summary(survfitdata, time=time)$upp)*100,2), 
             cuminc_upp = round((1 - summary(survfitdata, time=time)$low)*100,2), 
             scenario = scenario)
}

# risk by ages 18,25,30 for 2015
all_results_risk <- list()
years <- c(2015,2021) 
outcomes <- c("affek","epi")
for (year in years) {
  year_breaks <- if(year == 2015){c(18,25,30)} else {c(18,25,30,40)}
  for (out in outcomes){
    cat("Year:", year," , ", "Outcome:", out,  "\n")
    all_results_risk <- 
      append(all_results_risk, list(
        res_times(get(paste0("survfit_full_",out,year)),out,year,year_breaks, "full"),
        res_times(get(paste0("survfit_sc_",out,year)),out,year, year_breaks, "sc"),
        res_times(get(paste0("survfit_sc_",out,"_",out,year)),out,year, year_breaks, "sc_dis"),
        res_times(get(paste0("survfit_sc_any_",out,year)),out,year, year_breaks, "sc_any")
      ))
    cat("\n") # blank line
  }
}
comb_res_risk <- bind_rows(all_results_risk)

comb_res_risk %>% 
  filter(year == 2015) %>% 
  write.xlsx("risks_FU2015_aff_epi_test.xlsx",overwrite=TRUE)

comb_res_risk %>% 
  filter(year == 2021) %>% 
  write.xlsx("risks_FU2021_aff_epi_test.xlsx",overwrite=TRUE)

# Smoothed curves using GAM models predict at 500 values over age 
# 2015
# Affective disorder 
# Full cohort
gam_surv_full_affek2015 <- gam(p ~ s(time, k=10), data = survest_full_affek2015, method="REML")
gam_surv_ul_full_affek2015 <- gam(conf.high ~ s(time, k=10), data = survest_full_affek2015, method="REML")
gam_surv_ll_full_affek2015 <- gam(conf.low ~ s(time, k=10), data = survest_full_affek2015, method="REML")
# new time points for plotting
time2015 <- data.frame(time=seq(1, max(survfit_full_affek2015$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_full_affek2015 <- predict(gam_surv_full_affek2015, newdata = time2015, type = "response")
surv_p_smooth_ul_full_affek2015 <- predict(gam_surv_ul_full_affek2015, newdata = time2015, type = "response")
surv_p_smooth_ll_full_affek2015 <- predict(gam_surv_ll_full_affek2015, newdata = time2015, type = "response")
# SC
gam_surv_sc_affek2015 <- gam(p ~ s(time, k=10), data = survest_sc_affek2015, method="REML")
gam_surv_ul_sc_affek2015 <- gam(conf.high ~ s(time, k=10), data = survest_sc_affek2015, method="REML")
gam_surv_ll_sc_affek2015 <- gam(conf.low ~ s(time, k=10), data = survest_sc_affek2015, method="REML")
# new time points for plotting
time2015 <- data.frame(time=seq(1, max(survfit_sc_affek2015$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_affek2015 <- predict(gam_surv_sc_affek2015, newdata = time2015, type = "response")
surv_p_smooth_ul_sc_affek2015 <- predict(gam_surv_ul_sc_affek2015, newdata = time2015, type = "response")
surv_p_smooth_ll_sc_affek2015 <- predict(gam_surv_ll_sc_affek2015, newdata = time2015, type = "response")
# SC+affek 
gam_surv_sc_affek_affek2015 <- gam(p ~ s(time, k=10), data = survest_sc_affek_affek2015, method="REML")
gam_surv_ul_sc_affek_affek2015 <- gam(conf.high ~ s(time, k=10), data = survest_sc_affek_affek2015, method="REML")
gam_surv_ll_sc_affek_affek2015 <- gam(conf.low ~ s(time, k=10), data = survest_sc_affek_affek2015, method="REML")
# new time points for plotting
time2015 <- data.frame(time=seq(1, max(survfit_sc_affek_affek2015$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_affek_affek2015 <- predict(gam_surv_sc_affek_affek2015, newdata = time2015, type = "response")
surv_p_smooth_ul_sc_affek_affek2015 <- predict(gam_surv_ul_sc_affek_affek2015, newdata = time2015, type = "response")
surv_p_smooth_ll_sc_affek_affek2015 <- predict(gam_surv_ll_sc_affek_affek2015, newdata = time2015, type = "response")
# SC_any 
gam_surv_sc_any_affek2015 <- gam(p ~ s(time, k=10), data = survest_sc_any_affek2015, method="REML")
gam_surv_ul_sc_any_affek2015 <- gam(conf.high ~ s(time, k=10), data = survest_sc_any_affek2015, method="REML")
gam_surv_ll_sc_any_affek2015 <- gam(conf.low ~ s(time, k=10), data = survest_sc_any_affek2015, method="REML")
# new time points for plotting
time2015 <- data.frame(time=seq(1, max(survfit_sc_any_affek2015$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_any_affek2015 <- predict(gam_surv_sc_any_affek2015, newdata = time2015, type = "response")
surv_p_smooth_ul_sc_any_affek2015 <- predict(gam_surv_ul_sc_any_affek2015, newdata = time2015, type = "response")
surv_p_smooth_ll_sc_any_affek2015 <- predict(gam_surv_ll_sc_any_affek2015, newdata = time2015, type = "response")

# all predicted values for affek2015
pred_risk_affek2015_all <- as.data.frame(cbind(
  age = time2015$time,
  est_full = (1-p_smooth_full_affek2015)*100,
  ci_l_full=(1-surv_p_smooth_ul_full_affek2015)*100,
  ci_u_full=(1-surv_p_smooth_ll_full_affek2015)*100,
  est_sc = (1-p_smooth_sc_affek2015)*100,
  ci_l_sc=(1-surv_p_smooth_ul_sc_affek2015)*100,
  ci_u_sc=(1-surv_p_smooth_ll_sc_affek2015)*100,
  est_sc_affek = (1-p_smooth_sc_affek_affek2015)*100,
  ci_l_sc_affek=(1-surv_p_smooth_ul_sc_affek_affek2015)*100,
  ci_u_sc_affek=(1-surv_p_smooth_ll_sc_affek_affek2015)*100,
  est_sc_any = (1-p_smooth_sc_any_affek2015)*100,
  ci_l_sc_any=(1-surv_p_smooth_ul_sc_any_affek2015)*100,
  ci_u_sc_any=(1-surv_p_smooth_ll_sc_any_affek2015)*100
))
pred_risk_affek2015_all <- pred_risk_affek2015_all %>% mutate(year = "2015", disorder = "affek")

# Epilepsy
# Full cohort
gam_surv_full_epi2015 <- gam(p ~ s(time, k=10), data = survest_full_epi2015, method="REML")
gam_surv_ul_full_epi2015 <- gam(conf.high ~ s(time, k=10), data = survest_full_epi2015, method="REML")
gam_surv_ll_full_epi2015 <- gam(conf.low ~ s(time, k=10), data = survest_full_epi2015, method="REML")
# new time points for plotting
time2015 <- data.frame(time=seq(1, max(survfit_full_epi2015$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_full_epi2015 <- predict(gam_surv_full_epi2015, newdata = time2015, type = "response")
surv_p_smooth_ul_full_epi2015 <- predict(gam_surv_ul_full_epi2015, newdata = time2015, type = "response")
surv_p_smooth_ll_full_epi2015 <- predict(gam_surv_ll_full_epi2015, newdata = time2015, type = "response")

# SC
gam_surv_sc_epi2015 <- gam(p ~ s(time, k=10), data = survest_sc_epi2015, method="REML")
gam_surv_ul_sc_epi2015 <- gam(conf.high ~ s(time, k=10), data = survest_sc_epi2015, method="REML")
gam_surv_ll_sc_epi2015 <- gam(conf.low ~ s(time, k=10), data = survest_sc_epi2015, method="REML")
# new time points for plotting
time2015 <- data.frame(time=seq(1, max(survfit_sc_epi2015$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_epi2015 <- predict(gam_surv_sc_epi2015, newdata = time2015, type = "response")
surv_p_smooth_ul_sc_epi2015 <- predict(gam_surv_ul_sc_epi2015, newdata = time2015, type = "response")
surv_p_smooth_ll_sc_epi2015 <- predict(gam_surv_ll_sc_epi2015, newdata = time2015, type = "response")
# SC+epi
gam_surv_sc_epi_epi2015 <- gam(p ~ s(time, k=10), data = survest_sc_epi_epi2015, method="REML")
gam_surv_ul_sc_epi_epi2015 <- gam(conf.high ~ s(time, k=10), data = survest_sc_epi_epi2015, method="REML")
gam_surv_ll_sc_epi_epi2015 <- gam(conf.low ~ s(time, k=10), data = survest_sc_epi_epi2015, method="REML")
# new time points for plotting
time2015 <- data.frame(time=seq(1, max(survfit_sc_epi_epi2015$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_epi_epi2015 <- predict(gam_surv_sc_epi_epi2015, newdata = time2015, type = "response")
surv_p_smooth_ul_sc_epi_epi2015 <- predict(gam_surv_ul_sc_epi_epi2015, newdata = time2015, type = "response")
surv_p_smooth_ll_sc_epi_epi2015 <- predict(gam_surv_ll_sc_epi_epi2015, newdata = time2015, type = "response")
# SC_any
gam_surv_sc_any_epi2015 <- gam(p ~ s(time, k=10), data = survest_sc_any_epi2015, method="REML")
gam_surv_ul_sc_any_epi2015 <- gam(conf.high ~ s(time, k=10), data = survest_sc_any_epi2015, method="REML")
gam_surv_ll_sc_any_epi2015 <- gam(conf.low ~ s(time, k=10), data = survest_sc_any_epi2015, method="REML")
# new time points for plotting
time2015 <- data.frame(time=seq(1, max(survfit_sc_any_epi2015$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_any_epi2015 <- predict(gam_surv_sc_any_epi2015, newdata = time2015, type = "response")
surv_p_smooth_ul_sc_any_epi2015 <- predict(gam_surv_ul_sc_any_epi2015, newdata = time2015, type = "response")
surv_p_smooth_ll_sc_any_epi2015 <- predict(gam_surv_ll_sc_any_epi2015, newdata = time2015, type = "response")

# all predicted values for epi2015
pred_risk_epi2015_all <- as.data.frame(cbind(
  age = time2015$time,
  est_full = (1-p_smooth_full_epi2015)*100,
  ci_l_full=(1-surv_p_smooth_ul_full_epi2015)*100,
  ci_u_full=(1-surv_p_smooth_ll_full_epi2015)*100,
  est_sc = (1-p_smooth_sc_epi2015)*100,
  ci_l_sc=(1-surv_p_smooth_ul_sc_epi2015)*100,
  ci_u_sc=(1-surv_p_smooth_ll_sc_epi2015)*100,
  est_sc_epi = (1-p_smooth_sc_epi_epi2015)*100,
  ci_l_sc_epi=(1-surv_p_smooth_ul_sc_epi_epi2015)*100,
  ci_u_sc_epi=(1-surv_p_smooth_ll_sc_epi_epi2015)*100,
  est_sc_any = (1-p_smooth_sc_any_epi2015)*100,
  ci_l_sc_any=(1-surv_p_smooth_ul_sc_any_epi2015)*100,
  ci_u_sc_any=(1-surv_p_smooth_ll_sc_any_epi2015)*100
))

pred_risk_epi2015_all <- pred_risk_epi2015_all %>% mutate(year = "2015", disorder = "epi")

#2021
# Affective disorder
# Full cohort
gam_surv_full_affek2021 <- gam(p ~ s(time, k=10), data = survest_full_affek2021, method="REML")
gam_surv_ul_full_affek2021 <- gam(conf.high ~ s(time, k=10), data = survest_full_affek2021, method="REML")
gam_surv_ll_full_affek2021 <- gam(conf.low ~ s(time, k=10), data = survest_full_affek2021, method="REML")
# new time points for plotting
time2021 <- data.frame(time=seq(1, max(survfit_full_affek2021$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_full_affek2021 <- predict(gam_surv_full_affek2021, newdata = time2021, type = "response")
surv_p_smooth_ul_full_affek2021 <- predict(gam_surv_ul_full_affek2021, newdata = time2021, type = "response")
surv_p_smooth_ll_full_affek2021 <- predict(gam_surv_ll_full_affek2021, newdata = time2021, type = "response")

# SC
gam_surv_sc_affek2021 <- gam(p ~ s(time, k=10), data = survest_sc_affek2021, method="REML")
gam_surv_ul_sc_affek2021 <- gam(conf.high ~ s(time, k=10), data = survest_sc_affek2021, method="REML")
gam_surv_ll_sc_affek2021 <- gam(conf.low ~ s(time, k=10), data = survest_sc_affek2021, method="REML")
# new time points for plotting
time2021 <- data.frame(time=seq(1, max(survfit_sc_affek2021$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_affek2021 <- predict(gam_surv_sc_affek2021, newdata = time2021, type = "response")
surv_p_smooth_ul_sc_affek2021 <- predict(gam_surv_ul_sc_affek2021, newdata = time2021, type = "response")
surv_p_smooth_ll_sc_affek2021 <- predict(gam_surv_ll_sc_affek2021, newdata = time2021, type = "response")
# SC+affek
gam_surv_sc_affek_affek2021 <- gam(p ~ s(time, k=10), data = survest_sc_affek_affek2021, method="REML")
gam_surv_ul_sc_affek_affek2021 <- gam(conf.high ~ s(time, k=10), data = survest_sc_affek_affek2021, method="REML")
gam_surv_ll_sc_affek_affek2021 <- gam(conf.low ~ s(time, k=10), data = survest_sc_affek_affek2021, method="REML")
# new time points for plotting
time2021 <- data.frame(time=seq(1, max(survfit_sc_affek_affek2021$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_affek_affek2021 <- predict(gam_surv_sc_affek_affek2021, newdata = time2021, type = "response")
surv_p_smooth_ul_sc_affek_affek2021 <- predict(gam_surv_ul_sc_affek_affek2021, newdata = time2021, type = "response")
surv_p_smooth_ll_sc_affek_affek2021 <- predict(gam_surv_ll_sc_affek_affek2021, newdata = time2021, type = "response")
# SC_any
gam_surv_sc_any_affek2021 <- gam(p ~ s(time, k=10), data = survest_sc_any_affek2021, method="REML")
gam_surv_ul_sc_any_affek2021 <- gam(conf.high ~ s(time, k=10), data = survest_sc_any_affek2021, method="REML")
gam_surv_ll_sc_any_affek2021 <- gam(conf.low ~ s(time, k=10), data = survest_sc_any_affek2021, method="REML")
# new time points for plotting
time2021 <- data.frame(time=seq(1, max(survfit_sc_any_affek2021$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_any_affek2021 <- predict(gam_surv_sc_any_affek2021, newdata = time2021, type = "response")
surv_p_smooth_ul_sc_any_affek2021 <- predict(gam_surv_ul_sc_any_affek2021, newdata = time2021, type = "response")
surv_p_smooth_ll_sc_any_affek2021 <- predict(gam_surv_ll_sc_any_affek2021, newdata = time2021, type = "response")

# all predicted values for affek2021
pred_risk_affek2021_all <- as.data.frame(cbind(
  age = time2021$time,
  est_full = (1-p_smooth_full_affek2021)*100,
  ci_l_full=(1-surv_p_smooth_ul_full_affek2021)*100,
  ci_u_full=(1-surv_p_smooth_ll_full_affek2021)*100,
  est_sc = (1-p_smooth_sc_affek2021)*100,
  ci_l_sc=(1-surv_p_smooth_ul_sc_affek2021)*100,
  ci_u_sc=(1-surv_p_smooth_ll_sc_affek2021)*100,
  est_sc_affek = (1-p_smooth_sc_affek_affek2021)*100,
  ci_l_sc_affek=(1-surv_p_smooth_ul_sc_affek_affek2021)*100,
  ci_u_sc_affek=(1-surv_p_smooth_ll_sc_affek_affek2021)*100,
  est_sc_any = (1-p_smooth_sc_any_affek2021)*100,
  ci_l_sc_any=(1-surv_p_smooth_ul_sc_any_affek2021)*100,
  ci_u_sc_any=(1-surv_p_smooth_ll_sc_any_affek2021)*100
))

pred_risk_affek2021_all <- pred_risk_affek2021_all %>% mutate(year = "2021", disorder = "affek")

# Epilepsy
# Full cohort
gam_surv_full_epi2021 <- gam(p ~ s(time, k=10), data = survest_full_epi2021, method="REML")
gam_surv_ul_full_epi2021 <- gam(conf.high ~ s(time, k=10), data = survest_full_epi2021, method="REML")
gam_surv_ll_full_epi2021 <- gam(conf.low ~ s(time, k=10), data = survest_full_epi2021, method="REML")
# new time points for plotting
time2021 <- data.frame(time=seq(1, max(survfit_full_epi2021$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_full_epi2021 <- predict(gam_surv_full_epi2021, newdata = time2021, type = "response")
surv_p_smooth_ul_full_epi2021 <- predict(gam_surv_ul_full_epi2021, newdata = time2021, type = "response")
surv_p_smooth_ll_full_epi2021 <- predict(gam_surv_ll_full_epi2021, newdata = time2021, type = "response")

# SC
gam_surv_sc_epi2021 <- gam(p ~ s(time, k=10), data = survest_sc_epi2021, method="REML")
gam_surv_ul_sc_epi2021 <- gam(conf.high ~ s(time, k=10), data = survest_sc_epi2021, method="REML")
gam_surv_ll_sc_epi2021 <- gam(conf.low ~ s(time, k=10), data = survest_sc_epi2021, method="REML")
# new time points for plotting
time2021 <- data.frame(time=seq(1, max(survfit_sc_epi2021$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_epi2021 <- predict(gam_surv_sc_epi2021, newdata = time2021, type = "response")
surv_p_smooth_ul_sc_epi2021 <- predict(gam_surv_ul_sc_epi2021, newdata = time2021, type = "response")
surv_p_smooth_ll_sc_epi2021 <- predict(gam_surv_ll_sc_epi2021, newdata = time2021, type = "response")

# SC+epi
gam_surv_sc_epi_epi2021 <- gam(p ~ s(time, k=10), data = survest_sc_epi_epi2021, method="REML")
gam_surv_ul_sc_epi_epi2021 <- gam(conf.high ~ s(time, k=10), data = survest_sc_epi_epi2021, method="REML")
gam_surv_ll_sc_epi_epi2021 <- gam(conf.low ~ s(time, k=10), data = survest_sc_epi_epi2021, method="REML")
# new time points for plotting
time2021 <- data.frame(time=seq(1, max(survfit_sc_epi_epi2021$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_epi_epi2021 <- predict(gam_surv_sc_epi_epi2021, newdata = time2021, type = "response")
surv_p_smooth_ul_sc_epi_epi2021 <- predict(gam_surv_ul_sc_epi_epi2021, newdata = time2021, type = "response")
surv_p_smooth_ll_sc_epi_epi2021 <- predict(gam_surv_ll_sc_epi_epi2021, newdata = time2021, type = "response")
# SC_any
gam_surv_sc_any_epi2021 <- gam(p ~ s(time, k=10), data = survest_sc_any_epi2021, method="REML")
gam_surv_ul_sc_any_epi2021 <- gam(conf.high ~ s(time, k=10), data = survest_sc_any_epi2021, method="REML")
gam_surv_ll_sc_any_epi2021 <- gam(conf.low ~ s(time, k=10), data = survest_sc_any_epi2021, method="REML")
# new time points for plotting
time2021 <- data.frame(time=seq(1, max(survfit_sc_any_epi2021$time), length.out = 500))
# predict smoothed inverse survival prob and CI using GAM model
p_smooth_sc_any_epi2021 <- predict(gam_surv_sc_any_epi2021, newdata = time2021, type = "response")
surv_p_smooth_ul_sc_any_epi2021 <- predict(gam_surv_ul_sc_any_epi2021, newdata = time2021, type = "response")
surv_p_smooth_ll_sc_any_epi2021 <- predict(gam_surv_ll_sc_any_epi2021, newdata = time2021, type = "response")

# all predicted values for epi2021
pred_risk_epi2021_all <- as.data.frame(cbind(
  age = time2021$time,
  est_full = (1-p_smooth_full_epi2021)*100,
  ci_l_full=(1-surv_p_smooth_ul_full_epi2021)*100,
  ci_u_full=(1-surv_p_smooth_ll_full_epi2021)*100,
  est_sc = (1-p_smooth_sc_epi2021)*100,
  ci_l_sc=(1-surv_p_smooth_ul_sc_epi2021)*100,
  ci_u_sc=(1-surv_p_smooth_ll_sc_epi2021)*100,
  est_sc_epi = (1-p_smooth_sc_epi_epi2021)*100,
  ci_l_sc_epi=(1-surv_p_smooth_ul_sc_epi_epi2021)*100,
  ci_u_sc_epi=(1-surv_p_smooth_ll_sc_epi_epi2021)*100,
  est_sc_any = (1-p_smooth_sc_any_epi2021)*100,
  ci_l_sc_any=(1-surv_p_smooth_ul_sc_any_epi2021)*100,
  ci_u_sc_any=(1-surv_p_smooth_ll_sc_any_epi2021)*100
))

pred_risk_epi2021_all <- pred_risk_epi2021_all %>% mutate(year = "2021", disorder = "epi")

pred_risk_all <- bind_rows(pred_risk_affek2015_all,pred_risk_affek2021_all,
                           pred_risk_epi2015_all,pred_risk_epi2021_all)

write_csv(pred_risk_all,"pred_risk_aff_epi.csv")

# Figure 3 - plot cumulative risks across age for Affective, epilepsy ----------

# read dataframe including gam-smoothed cum inc.
pred_risk_aff_epi <- read_csv("pred_risk_aff_epi.csv")

# one variable with estimate for scenario sc_disorder
pred_risk_aff_epi <- pred_risk_aff_epi %>%
  mutate(
    est_sc_dis = case_when(
      disorder == "affek" ~ est_sc_affek,
      disorder == "epi" ~ est_sc_epi),
    ci_l_sc_dis = case_when(disorder == "affek" ~ ci_l_sc_affek,
                            disorder == "epi" ~ ci_l_sc_epi),
    ci_u_sc_dis = case_when(disorder == "affek" ~ ci_u_sc_affek,
                            disorder == "epi" ~ ci_u_sc_epi)
  )

# plot risk --------------------------------------------------------------------
plot3 <- pred_risk_aff_epi %>%
  mutate(disorder = factor(disorder,
                           levels = c("affek","epi"), 
                           labels = c("Affective disorder","Epilepsy")
  )) %>%
  mutate(year = factor(year,
                       levels = c("2015","2021"), 
                       labels = c("Follow-up by 2015","Follow-up by 2021")
  )) %>%
  ggplot( aes(x=age)) +
  geom_line(aes(y=est_full, color = "Full cohort"), lwd=0.3) +
  geom_ribbon(aes(ymin = ci_l_full, ymax = ci_u_full, fill="Full cohort"), alpha = 0.6, show.legend=FALSE) +
  geom_line(aes(y=est_sc, color = "Subcohort"), lwd=0.3) +
  geom_ribbon(aes(ymin = ci_l_sc, ymax = ci_u_sc, fill="Subcohort"), alpha = 0.1, show.legend=FALSE) +
  geom_line(aes(y=est_sc_dis, color = "Subcohort and outcome"), lwd=0.3) +
  geom_ribbon(aes(ymin = ci_l_sc_dis, ymax = ci_u_sc_dis, fill="Subcohort and outcome"), alpha = 0.1, show.legend=FALSE) +
  geom_line(aes(y=est_sc_any, color = "Subcohort and all cases"), lwd=0.3) +
  geom_ribbon(aes(ymin = ci_l_sc_any, ymax = ci_u_sc_any, fill="Subcohort and all cases"), alpha = 0.1, show.legend=FALSE) +
  #NOTE: Never version of the facet_grid can have: axes="all", axis.labels = "marginal" 
  # but for now we need to use the lemon package facet_rep_grid to get that effect.
  lemon::facet_rep_grid(disorder ~ year,
                        , scale = "free", switch = "y",space = "free_x") +#,axes="all"
  theme_minimal() +
  scale_y_continuous(sec.axis = dup_axis(labels = NULL)) + 
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle=0, size=8),
        strip.text.x.top = element_text(size=10),
        legend.position = "bottom",
        legend.text = element_text(size=8),
        axis.title.y.left = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.y.left = element_line(),
        axis.text.y.left = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.line=element_line(color="black"),
        axis.line.y.right = element_blank()
  ) +
  scale_x_continuous(breaks = seq(0,40,by=5)) +
  #                    ) + 
  labs(x="Age (years)", y="Risk (%)") +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks.x = element_line(),
        axis.title.y.right = element_text(size=8),
        axis.title.x = element_text(size=8),
  ) +
  scale_color_manual(name="", values=c("Full cohort" = "black","Subcohort"= "red", "Subcohort and outcome"="chartreuse2","Subcohort and all cases"="blue"),
                     breaks = c("Full cohort","Subcohort","Subcohort and outcome","Subcohort and all cases")) +
  scale_fill_manual(values = c("Full cohort" = "grey", "Subcohort" = "red", "Subcohort and outcome"="chartreuse2", "Subcohort and all cases"="blue"))

plot3

ggsave("risk_aff_epi.svg", plot = plot3, width = 8, height = 4)


