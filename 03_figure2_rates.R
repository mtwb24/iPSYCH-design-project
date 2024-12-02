# setwd("E:/Data/workdata/708279/WOI/iPSYCH2015design/work_temp")
# ==============================================================================
# Title: Smoothed rates across age (Figure 2)
# Author: TW
# Date: Nov 2024
# Description: Estimating GAM poisson rates and plot across age for affective
# disorder and epilepsy
# ==============================================================================

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(survival)
library(Epi)
library(mgcv) 
library(ggplot2)
library(lemon) 

# Import data
setwd("...")
studybase <- read_csv("studybase.csv")

## time to event - function ----------------------------------------------------
# select variables, define survival data, indicator and date of event/censoring
studybase_tte <- function(data,eventdvar, age_tstart, lastfu, emidate, sample){
  data %>%
    select(xnr,fdato.x,{{eventdvar}}, death_d,{{emidate}}, born1981_2005, born2006_2008, {{lastfu}},kontrol2015i, affek2015i, 
           case_orig, case_any, full, w_full, sc_any, w_sc_any, paste0(sample), paste0("w_",sample), sc, w_sc, kqn.x) %>%  
    mutate(tstop = pmin({{eventdvar}}, {{lastfu}}, death_d, {{emidate}}, 
                        na.rm = TRUE)) %>%
    mutate(event = ifelse({{eventdvar}}>tstop | is.na({{eventdvar}}),0,1),
           age_tstart = age_tstart,
           age_tstop = (unclass(tstop) - unclass(fdato.x) ) /365.24,
           id = xnr,
           death = death_d == tstop
           ) 
}

## time to event - data frames using function ----------------------------------
# Note that date of first affective disorder diagnosis after age 
# 10 and is given in the variable affek2015 (iPSYCH design variable) for 
# diagnoses by 2015 and in the variable aff10_1d (updated information from 
# latest version of the Danish National Patient Register) for diagnoses by 2021. 
# Note: individuals are followed from age 1 in all scenarios.

# affective disorder, follow-up by 2015 and 2021
studybase_aff2015 <- studybase_tte(studybase,affek2015,age_tstart=1,last2015,emi10y_d,sample="sc_affek")
studybase_aff2021 <- studybase_tte(studybase,aff10_1d, age_tstart=1,last2021,emi10y_d,sample="sc_affek")
# epilepsy, follow-up by 2015 and 2021
studybase_epi2015 <- studybase_tte(studybase,epi_1d,age_tstart=1,last2015,emi1y_d,sample="sc_epi")
studybase_epi2021 <- studybase_tte(studybase,epi_1d, age_tstart=1,last2021,emi1y_d,sample="sc_epi")

## agesplit - function ---------------------------------------------------------
# function surv_agesplit: creating survival object and split data by age groups,
# default is cut age from 1 to 40 in steps of 1
surv_agesplit <- function(data,event,tstart, tstop,age_breaks = c(seq(1,40, by=1))){
  survSplit(Surv(age_tstart, age_tstop, event) ~ ., data = data, cut = age_breaks) 
}
## agesplit - data frames using function ---------------------------------------
# affective disorder, follow-up by 2015 and 2021
studybase_agesplit_aff2015 <- surv_agesplit(studybase_aff2015, 'event', 'age_tstart', 'age_tstop', age_breaks = c(seq(1,34, by=1)))
studybase_agesplit_aff2021 <- surv_agesplit(studybase_aff2021, 'event', 'age_tstart', 'age_tstop', age_breaks = c(seq(1,40, by=1)))
# epilepsy, follow-up by 2015 and 2021
studybase_agesplit_epi2015 <- surv_agesplit(studybase_epi2015, 'event', 'age_tstart', 'age_tstop', age_breaks = c(seq(1,34, by=1)))
studybase_agesplit_epi2021 <- surv_agesplit(studybase_epi2021, 'event', 'age_tstart', 'age_tstop', age_breaks = c(seq(1,40, by=1)))

# GAM (generalized additive models) Poisson regression rates by age by 2015 ----
# prediction data frame
age_pred2015 <- data.frame(age_tstop=seq(1, max(studybase_aff2015$age_tstop), length.out=500))
# affective disorder by 2015
# full cohort - for smoothed incidence rates curves, Carsten Bendixen notation, s. 135
pois_gam_full_aff2015 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                             family = poisreg, data = studybase_agesplit_aff2015)
pred_gam_full_aff2015 <- ci.pred(pois_gam_full_aff2015,newdata=age_pred2015) # predict fitted values over age (Wald CIs)
matshade(age_pred2015$age_tstop,pred_gam_full_aff2015,plot=TRUE,lwd=2, xlab="Age", ylab="Incidence rate/1000 person-years")
# subcohort
pois_gam_sc_affek2015 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                             family = poisreg, data = studybase_agesplit_aff2015, subset=sc==1)
pred_gam_sc_affek2015 <- ci.pred(pois_gam_sc_affek2015,newdata=age_pred2015)
# subcohort + affec, OBS: not robust confidence intervals, plotting est only
pois_gam_sc_affek_aff2015 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                                 family = poisreg, data = studybase_agesplit_aff2015, subset=sc_affek==1, weights=w_sc_affek)
pred_gam_sc_affek_aff2015 <- ci.pred(pois_gam_sc_affek_aff2015,newdata=age_pred2015)
# subcohort + any
pois_gam_sc_any_aff2015 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                               family = poisreg, data = studybase_agesplit_aff2015, subset=sc_any==1, weights=w_sc_any)
pred_gam_sc_any_aff2015 <- ci.pred(pois_gam_sc_any_aff2015,newdata=age_pred2015)
# dataframe with all estimates, no 95% CIs for case-cohort samples
pred_gam_aff2015_all <-
  as.data.frame(
    cbind(
      age = age_pred2015$age_tstop,
      rate_full = pred_gam_full_aff2015[, "Estimate"],
      ci_l_full = pred_gam_full_aff2015[, "2.5%"],
      ci_u_full = pred_gam_full_aff2015[, "97.5%"],
      rate_sc = pred_gam_sc_affek2015[, "Estimate"],
      ci_l_sc = pred_gam_sc_affek2015[, "2.5%"],
      ci_u_sc = pred_gam_sc_affek2015[, "97.5%"],
      rate_sc_affek = pred_gam_sc_affek_aff2015[, "Estimate"],
      rate_sc_any = pred_gam_sc_any_aff2015[, "Estimate"]
    )
  )
pred_gam_aff2015_all <- pred_gam_aff2015_all %>%
  mutate(year = "2015", disorder = "aff")

# epilepsy by 2015 
pois_gam_full_epi2015 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                             family = poisreg, data = studybase_agesplit_epi2015)
pred_gam_full_epi2015 <- ci.pred(pois_gam_full_epi2015,newdata=age_pred2015)
pois_gam_sc_epi2015 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                           family = poisreg, data = studybase_agesplit_epi2015, subset=sc==1)
pred_gam_sc_epi2015 <- ci.pred(pois_gam_sc_epi2015,newdata=age_pred2015)
pois_gam_sc_epi_epi2015 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                               family = poisreg, data = studybase_agesplit_epi2015, subset=sc_epi==1, weights=w_sc_epi)
pred_gam_sc_epi_epi2015 <- ci.pred(pois_gam_sc_epi_epi2015,newdata=age_pred2015)
pois_gam_sc_any_epi2015 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                               family = poisreg, data = studybase_agesplit_epi2015, subset=sc_any==1, weights=w_sc_any)
pred_gam_sc_any_epi2015 <- ci.pred(pois_gam_sc_any_epi2015,newdata=age_pred2015)
# dataframe with all estimates, no 95% CIs for case-cohort samples
pred_gam_epi2015_all <- as.data.frame(
  cbind(
    age = age_pred2015$age_tstop,
    rate_full = pred_gam_full_epi2015[, "Estimate"],
    ci_l_full = pred_gam_full_epi2015[, "2.5%"],
    ci_u_full = pred_gam_full_epi2015[, "97.5%"],
    rate_sc = pred_gam_sc_epi2015[, "Estimate"],
    ci_l_sc = pred_gam_sc_epi2015[, "2.5%"],
    ci_u_sc = pred_gam_sc_epi2015[, "97.5%"],
    rate_sc_epi = pred_gam_sc_epi_epi2015[, "Estimate"],
    rate_sc_any = pred_gam_sc_any_epi2015[, "Estimate"]
  )
)
pred_gam_epi2015_all <- pred_gam_epi2015_all %>%
  mutate(year = "2015", disorder = "epi")

# GAM (generalized additive models) Poisson regression rates by age by 2021 ----
# prediction data frame, consider change to length 100 or 0.25 steps.
age_pred2021 <- data.frame(age_tstop=seq(1, max(studybase_aff2021$age_tstop), length.out=500))

# GAM rates, Affective disorder by 2021 (Full cohort, SC, sc_affek, all) (notice: warning when running full cohort)
pois_gam_aff2021 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                        family = poisreg, data = studybase_agesplit_aff2021)
pred_gam_full_aff2021 <- ci.pred(pois_gam_aff2021,newdata=age_pred2021)
matshade(age_pred2021$age_tstop,pred_gam_full_aff2021,plot=TRUE,lwd=2, xlab="Age", ylab="Incidence rate/1000 person-years")
pois_gam_sc_affek2021 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                             family = poisreg, data = studybase_agesplit_aff2021, subset=sc==1)
pred_gam_sc_affek2021 <- ci.pred(pois_gam_sc_affek2021,newdata=age_pred2021)
matshade(age_pred2021$age_tstop,pred_gam_sc_affek2021,plot=TRUE,lwd=2)
pois_gam_sc_affek_aff2021 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                                 family = poisreg, data = studybase_agesplit_aff2021, subset=sc_affek==1, weights=w_sc_affek)
pred_gam_sc_affek_aff2021 <- ci.pred(pois_gam_sc_affek_aff2021,newdata=age_pred2021)
matshade(age_pred2021$age_tstop,pred_gam_sc_affek_aff2021,plot=TRUE,lwd=2, xlab="Age", ylab="Incidence rate/1000 person-years")
pois_gam_sc_any_aff2021 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                               family = poisreg, data = studybase_agesplit_aff2021, subset=sc_any==1, weights=w_sc_any)
pred_gam_sc_any_aff2021 <- ci.pred(pois_gam_sc_any_aff2021,newdata=age_pred2021)
matshade(age_pred2021$age_tstop,pred_gam_sc_any_aff2021,plot=TRUE,lwd=2, xlab="Age", ylab="Incidence rate/1000 person-years")
# dataframe with all estimates, no 95% CIs for case-cohort samples
pred_gam_aff2021_all <- as.data.frame(
  cbind(
    age = age_pred2021$age_tstop,
    rate_full = pred_gam_full_aff2021[, "Estimate"],
    ci_l_full = pred_gam_full_aff2021[, "2.5%"],
    ci_u_full = pred_gam_full_aff2021[, "97.5%"],
    rate_sc = pred_gam_sc_affek2021[, "Estimate"],
    ci_l_sc = pred_gam_sc_affek2021[, "2.5%"],
    ci_u_sc = pred_gam_sc_affek2021[, "97.5%"],
    rate_sc_affek = pred_gam_sc_affek_aff2021[, "Estimate"],
    rate_sc_any = pred_gam_sc_any_aff2021[, "Estimate"]
  )
)
pred_gam_aff2021_all <- pred_gam_aff2021_all %>% mutate(year = "2021", disorder = "aff")


# GAM rates, Epilepsy by 2021 (Full cohort, SC, sc_epi, all)
pois_gam_epi2021 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                        family = poisreg, data = studybase_agesplit_epi2021)
pred_gam_full_epi2021 <- ci.pred(pois_gam_epi2021,newdata=age_pred2021)
pois_gam_sc_epi2021 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                           family = poisreg, data = studybase_agesplit_epi2021, subset=sc==1)
pred_gam_sc_epi2021 <- ci.pred(pois_gam_sc_epi2021,newdata=age_pred2021)
pois_gam_sc_epi_epi2021 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                               family = poisreg, data = studybase_agesplit_epi2021, subset=sc_epi==1, weights=w_sc_epi)
pred_gam_sc_epi_epi2021 <- ci.pred(pois_gam_sc_epi_epi2021,newdata=age_pred2021)
pois_gam_sc_any_epi2021 <- gam(cbind(event,(age_tstop-age_tstart)/1000) ~ s(age_tstop,k=7), 
                               family = poisreg, data = studybase_agesplit_epi2021, subset=sc_any==1, weights=w_sc_any)
pred_gam_sc_any_epi2021 <- ci.pred(pois_gam_sc_any_epi2021,newdata=age_pred2021)

# dataframe with all estimates, no 95% CIs for case-cohort samples
pred_gam_epi2021_all <- as.data.frame(
  cbind(
    age = age_pred2021$age_tstop,
    rate_full = pred_gam_full_epi2021[, "Estimate"],
    ci_l_full = pred_gam_full_epi2021[, "2.5%"],
    ci_u_full = pred_gam_full_epi2021[, "97.5%"],
    rate_sc = pred_gam_sc_epi2021[, "Estimate"],
    ci_l_sc = pred_gam_sc_epi2021[, "2.5%"],
    ci_u_sc = pred_gam_sc_epi2021[, "97.5%"],
    rate_sc_epi = pred_gam_sc_epi_epi2021[, "Estimate"],
    rate_sc_any = pred_gam_sc_any_epi2021[, "Estimate"]
  )
)
pred_gam_epi2021_all <- pred_gam_epi2021_all %>% mutate(year = "2021", disorder = "epi")


pred_gam_all <- bind_rows(pred_gam_aff2015_all,pred_gam_aff2021_all,
                          pred_gam_epi2015_all,pred_gam_epi2021_all)

# combine disorder-specific rates into one variable
pred_gam_all <- pred_gam_all %>%
  mutate(rate_sc_dis = case_when(
    disorder=="aff" ~ rate_sc_affek,
    disorder=="epi" ~ rate_sc_epi ))

## save GAM rates for plotting -------------------------------------------------
write_csv(pred_gam_all,"pred_gam_aff_epi.csv")

## Plot rates across age (Figure 2) --------------------------------------------

# plot rates across age for Affective, epilepsy, different samples and follow-up

# read dataframe including gam-smoothed poisson rates
pred_gam_all <- read_csv("pred_gam_aff_epi.csv")

plot2 <- pred_gam_all %>% filter(disorder == "aff" |
                                   disorder == "epi") %>%
  mutate(disorder = factor(
    disorder,
    levels = c("aff", "epi"),
    labels = c("Affective disorder", "Epilepsy")
  )) %>%
  mutate(year = factor(
    year,
    levels = c("2015", "2021"),
    labels = c("Follow-up by 2015", "Follow-up by 2021")
  )) %>%
  ggplot(aes(x = age)) +
  geom_line(aes(y = rate_full, color = "Full cohort"), lwd = 0.3) +
  geom_line(aes(y = rate_sc, color = "Subcohort"), lwd = 0.3) +
  geom_line(aes(y = rate_sc_dis, color = "Subcohort and outcome"), lwd =
              0.3) +
  geom_line(aes(y = rate_sc_any, color = "Subcohort and all cases"), lwd =
              0.3) +
  #NOTE: Never version of the facet_grid can have: axes="all", axis.labels = "marginal"
  # but for now we need to use the lemon package facet_rep_grid to get that effect.
  lemon::facet_rep_grid(disorder ~ year, # labeller = as_labeller(c(ans = )),
                        scale = "free",
                        switch = "y",
                        space = "free_x") + #,axes="all"
  theme_minimal() +
  scale_y_continuous(sec.axis = dup_axis(labels = NULL)) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 8),
    strip.text.x.top = element_text(size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    axis.title.y.left = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.ticks.y.left = element_line(),
    axis.text.y.left = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank()
  ) +
  scale_x_continuous(breaks = seq(0, 40, by = 5)) +
  labs(x = "Age (years)", y = "Incidence rate/1000 person-years") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(),
    axis.title.y.right = element_text(size = 8),
    axis.title.x = element_text(size = 8),
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "Full cohort" = "black",
      "Subcohort" = "red",
      "Subcohort and outcome" = "chartreuse2",
      "Subcohort and all cases" = "blue"
    ),
    breaks = c(
      "Full cohort",
      "Subcohort",
      "Subcohort and outcome",
      "Subcohort and all cases"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Full cohort" = "grey",
      "Subcohort" = "red",
      "Subcohort and outcome" = "chartreuse2",
      "Subcohort and all cases" = "blue"
    )
  )
plot2

ggsave(file.path("rates_aff_epi.svg"), plot = plot2, width = 8, height = 4, dpi = 300)

