# setwd("E:/Data/workdata/708279/WOI/iPSYCH2015design/work_temp")
# ==============================================================================
# Title: Hazard ratios using autism as time-dependent exposure (Figure 4)
# Author: TW
# Date: Nov 2024
# Description: Define time-varying exposure, estimate Hazard ratios from Cox 
# regression models with affective/epilepsy as outcomes and autism as exposure
# across four samples and for follow-up by 2015 and 2021.
# ==============================================================================

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(survival)
library(grid)
library(ggpubr)
library(patchwork)

# Import data
#setwd("...")
studybase <- read_csv("studybase.csv")

## time to event - function ----------------------------------------------------
# select variables, define survival data, indicator and date of event/censoring
studybase_tte <- function(data,eventdvar, age_tstart, lastfu, emidate, sample){
  data %>%
    select(xnr,fdato.x,{{eventdvar}}, death_d,{{emidate}}, born1981_2005, born2006_2008, {{lastfu}},kontrol2015i, affek2015i, 
           case_any, full, w_full, sc_any, w_sc_any, paste0(sample), paste0("w_",sample), sc, w_sc, kqn.x) %>%  
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
studybase_affek2015 <- studybase_tte(studybase,affek2015,age_tstart=1,last2015,emi10y_d,sample="sc_affek")
studybase_affek2021 <- studybase_tte(studybase,aff10_1d, age_tstart=1,last2021,emi10y_d,sample="sc_affek")
# epilepsy, follow-up by 2015 and 2021
studybase_epi2015 <- studybase_tte(studybase,epi_1d,age_tstart=1,last2015,emi1y_d,sample="sc_epi")
studybase_epi2021 <- studybase_tte(studybase,epi_1d, age_tstart=1,last2021,emi1y_d,sample="sc_epi")

# Time-dependent exposure: ASD
# followed from age 1y

# Cox models, time-dependent exposure: asd
# ASD as time-varying exposure
ASD <- studybase %>%
  select(xnr,fdato.x,asd_1d) %>%
  filter(!is.na(asd_1d))
# ASD times
asd_times <- ASD %>%
  mutate(asd_time = as.numeric((asd_1d - fdato.x)/365.24),
         asd_exp = 1) %>%
  select(xnr, asd_time, asd_exp)

# by outcome and year
all_results_asd <- list()

years <- c(2015,2021) 
outcomes <- c("affek","epi")
models <- list()
for (year in years) {
  for (out in outcomes){
    cat("Year:", year,",", "Exposure: asd,", "Outcome:", out,  "\n")
    # base data
    assign(paste0("base1_",out,year),get(paste0("studybase_",out,year)) %>% 
             select(xnr, event, age_tstart, age_tstop, fdato.x, kqn.x, 
                    full, sc, paste0("sc_",out), sc_any, paste0("w_sc_",out), w_sc_any) %>% 
             rename(tstart = age_tstart) %>%
             mutate(
               surv_time = age_tstop,
               sex = ifelse(kqn.x=="M",1,0)
             ))
    # base data with tstop need to be set up using tmerge
    base1 <- get(paste0("base1_",out,year))
    
    base2 <- tmerge(base1, base1, id=xnr, tstart = 1, tstop = age_tstop, event=event(surv_time,event))
    assign(paste0("base2_",out,year),base2)
    
    # merge base data with ADHD times
    df_asd <- tmerge(data1 = get(paste0("base2_",out,year)),
                     data2 = asd_times,
                     id = xnr,
                     exposure = tdc(asd_time,init=1)) 
    
    assign(paste0("df_asd_",out),df_asd)
    
    # Cox model
    m1 <- coxph(Surv(tstart,tstop,event) ~ exposure, data = get(paste0("df_asd_",out))) 
    m2 <- coxph(Surv(tstart,tstop,event) ~ exposure, data = get(paste0("df_asd_",out)), subset=sc==1)
    m3 <- coxph(Surv(tstart,tstop,event) ~ exposure + cluster(xnr), data = get(paste0("df_asd_",out)), subset=get(paste0("sc_",out))==1, weights = get(paste0("w_sc_",out)))
    m4 <- coxph(Surv(tstart,tstop,event) ~ exposure + cluster(xnr), data = get(paste0("df_asd_",out)), subset=sc_any==1, weights = w_sc_any)
    models[[as.character(year)]] <- list(m1=m1, m2=m2, m3=m3, m4=m4)
    s1 <- summary(m1)
    s2 <- summary(m2)
    s3 <- summary(m3)
    s4 <- summary(m4)
    # results
    res_asd <- list(Outcome = rep(out, times=4),
                    Year = rep(year, times=4),
                    Sample = c("Full cohort", "Subcohort", paste0("SC+",out),"SC+any"),
                    HR = c(s1$coefficients[1,2], s2$coefficients[1,2], s3$coefficients[1,2], s4$coefficients[1,2]),
                    ci_low = c(s1$conf.int[1,3], s2$conf.int[1,3], s3$conf.int[1,3], s4$conf.int[1,3]),
                    ci_upp = c(s1$conf.int[1,4], s2$conf.int[1,4], s3$conf.int[1,4], s4$conf.int[1,4]),
                    se = c(s1$coef[1,3], s2$coef[1,3], s3$coef[1,3], s4$coef[1,3]),
                    se_robust = c(s1$coef[1,3], s2$coef[1,3], s3$coef[1,4], s4$coef[1,4])) 
    
    all_results_asd <- append(all_results_asd, list(res_asd))
    cat("\n") # blank line
  }
}

HR_asd <- bind_rows(all_results_asd) 

# sample and group label variables
HR_asd <- HR_asd %>% mutate(
  sample = if_else(Sample %in% c("SC+affek","SC+epi"),"sc_disorder",Sample),
  sample = recode(sample, "Full cohort" = "full", "Subcohort" = "sc", "SC+any" = "sc_any"),
  grp = paste(Outcome,sample,sep="_"),
  # convert grp, sample, and Outcome variables to a factor with desired order (to obtain correct order in plot)
  grp = factor(grp,levels=unique(grp)),
  sample = factor(sample,levels=unique(sample)),
  disorder = factor(Outcome,levels=unique(Outcome)),
  # round estimates to (and display) 2 decimals, 3 decimals for standard errors
  est_lab = paste(sprintf("%.2f",HR), " (",sprintf("%.2f",ci_low),"-",sprintf("%.2f",ci_upp),")",sep=""),
  #se_robust <- as.numeric(se_robust),
  se_robust = paste(sprintf("%.3f",as.numeric(se_robust)))
)

# shorter grp var labels
HR_asd$grp2 <- HR_asd$grp
levels(HR_asd$grp2) <- 1:8
str(HR_asd$grp2)
# disorder label
HR_asd$disorder_lab <- HR_asd$disorder
levels(HR_asd$disorder_lab) <- c("Affective disorder","Epilepsy")

# subgroup df by year
HR_asd_2015 <- HR_asd %>% filter(Year=="2015") %>% 
  select(disorder, sample, grp, grp2, HR, ci_low, ci_upp, disorder_lab, est_lab, se_robust)
HR_asd_2021 <- HR_asd %>% filter(Year=="2021") %>% 
  select(disorder, sample, grp, grp2, HR, ci_low, ci_upp, disorder_lab, est_lab, se_robust)

# 1. create plot with estimates, using facet_wrap (by sample) for 2015 and 2021, respectively
forest_plot <- \(ds) {
  ggplot(data=ds, aes(x=HR,y=fct_rev(grp2),xmin=ci_low,xmax=ci_upp,col=sample)) +
    geom_vline(xintercept=1,linetype=2) +
    geom_pointrange(aes(), size=0.1)+
    geom_errorbar(width=0.2)+
    facet_wrap(~disorder,strip.position="left",ncol=1, ,scales="free_y")+
    ylab("")+
    xlab("HR") +
    labs(title = "") +
    theme_void()+
    theme(
      axis.ticks.y=element_blank(),
      axis.line.x=element_line(color="black"),
      axis.ticks= element_line(color="black"),
      axis.ticks.length = unit(0.1,"cm"),
      axis.title.x= element_text(size=8), # HR
      axis.text.x=element_text(size=8),  
      strip.text.y=element_blank(), 
      legend.position="none",
      legend.text = element_text(size = 8)
    )+
    scale_color_manual(guide = guide_legend(title = "Sample",reverse=FALSE), 
                       values = c("full"="black",
                                  "sc"="red",
                                  "sc_disorder"="chartreuse2",
                                  "sc_any"="blue"),
                       labels = c("full"="Full cohort",
                                  "sc"="Subcohort",
                                  "sc_disorder"="Subcohort and outcome", 
                                  "sc_any"="Subcohort and all cases")) +
    scale_x_continuous(breaks = c(0.5,1,seq(2,8, by=2)), 
                       limits=c(0.5,9), 
                       labels=c(0.5,1,
                       seq(2,8, by=2)), 
                       trans="log",
                       expand=c(0,0))
}

p2015 <- forest_plot(ds=HR_asd_2015)
p2015
p2021 <- forest_plot(ds=HR_asd_2021)
p2021

# disorder and estimate labs without headers (grobs 1 + 3)
p_text_groups <- HR_asd_2015 %>% 
  ggplot(aes(x = ""   ,y = disorder   , label = disorder_lab    )) +
  geom_text(angle=90,hjust = 0.5,vjust=4,size=4
  ) +
  facet_wrap(~disorder_lab,ncol=1, strip.position="left",scales="free_y")+
  scale_x_discrete(position = "top", labels = " ") +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    axis.text.x = element_text(hjust = 0),
    strip.text.y=element_blank()
  ) 
p_text_groups 


table_text <- \(ds){
  mutate(ds,se_robust = paste(se_robust)) %>%
    pivot_longer(cols = c("est_lab","se_robust"), names_to = "name_of_estimate") %>% 
    ggplot(aes(x = name_of_estimate,y = fct_rev(grp2), label = value )) +
    geom_text(hjust = 0,size=4) +
    facet_wrap(~disorder_lab, strip.position="left",ncol=1,scales="free_y"
    )+
    scale_x_discrete(position = "top", labels = c("HR (95% CI)","SE")) + 
    theme_void() +
    theme(
      axis.text.x = element_text(hjust = 0),
      strip.text.y=element_blank()
    )
} 
p_text2015 <- table_text(HR_asd_2015)
p_text2021 <- table_text(HR_asd_2021) 

str(p_text2015)

#Legend
leg <- ggpubr::get_legend(p2015+ theme(legend.title = element_blank()),position = "bottom")
ggleg <- ggpubr::as_ggplot(leg)

#Headline
title2015 <- textGrob("Follow-up until 2015", gp=gpar(fontsize=12))
title2021 <- textGrob("Follow-up until 2021", gp=gpar(fontsize=12))

# layout specification - position of plots notation: 
#1st row: blank slot (#), plots F and g each spanning two columns
#2nd row: five plots A-E
#3rd row: blank slot(#), H spanning four columns, centered in the row
layout <- "
#FFGG
ABCDE
#HHHH"

plot4 <- p_text_groups + p_text2015 + p2015 + p_text2021 + p2021  +
  title2015 + title2021 + 
  ggleg +
  plot_layout(design = layout,heights = c(0.1,1,0.2),widths = c(0.2,3,2,3,2) %>% (\(x) x/sum(x)) )

plot4

ggsave("HR_aff_epi_asd.svg", plot = plot4, width = 10, height = 4)