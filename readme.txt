# iPSYCH design paper

Github repository for code used to produce main figures and tables for the paper: "Population-Representative Inference for Primary and Secondary Outcomes in Extended Case-Cohort Designs: an Application of the iPSYCH Study"

In this study, we identify the entire source population for the iPSYCH sampl (the full cohort), nested subcohorts and case-cohort samples. We define average- and birth-cohort-specific inverse probability of sampling weights, and estimate age-specific rates, risks, and incidence rate ratios using crude poisson, GAM, and Cox regression modeling. We here use affective disorder as an example of a primary outcome, and epilepsy as an example of a secondary outcome. Estimation is done in the different samples (full cohort, subcohort, subcohort + affective disorder/epilepsy, and subcohort + all case groups = entire iPSYCH sample) and for different length of follow-up (end of 2015 and 2021, respectively).

### R-version 
R-4.3.2

### Libraries used

| Library       | Used in                         |
| ------------- | ------------------------------- |
| haven         | data management                 |
| tidiverse     | data management                 |
| openxlsx      | read and write Excel files      |
| survival      | data management and analysis    |
| Epi           | gam modeling                    |
| mgcv          | gam modeling                    |
| ggplot2       | creating figures                |
| lemon         | figure layout                   |
| broom         | data management                 |
| grid          | figure layout                   |
| ggpubr        | figure layout                   |
| LexisPlotR    | figure layout                   |

### Scripts overview

#### 01_data_management.R
This script defines samples and weights for primary and
secondary outcomes until 2015 and 2021 based on datasets including information
on the iPSYCH design and diagnoses from the National Patient Register.

Output file: studybase.csv

#### 02_table1_descriptives.R
This file calculates crude and upweighted numbers for different case-cohort samples in table 1.

Output file: upweight_tab1.xlsx

#### 03_figure2_rates.R
Estimating GAM poisson rates and plot across age for affective disorder and epilepsy (Figure 2).

Output file: rates_aff_epi.svg 

#### 04_figure3_risks.R
Cumulative incidence estimates (1-KM) across age for affective disorder (primary outcome) and epilepsy (secondary outcome) using GAM (generalized additive models) for smoothing. Across four samples and with follow-up by calendar years 2015 and 2021. Estimates calculated for ages 18, 25, 30, (40), and plotted over age (Figure 3).

Output files: risks_FU2015_aff_epi.xlsx, risks_FU2021_aff_epi.xlsx, risk_aff_epi.svg

#### 05_figure4_IRR.R
This file defines time-dependent exposure and estimate Hazard ratios from Cox 
regression models with affective/epilepsy as outcomes and autism as exposure (Figure 4).
Across four samples and for follow-up by 2015 and 2021.

Output file: HR_aff_epi_asd.svg

#### 06_sfigure1_lexis.R
Lexisplot illustrating iPSYCH cohorts and follow-up (sFigure 1)

Output file: sfig1_lexis_iPSYCH.png

### Availability of data and materials

Data for this study is property of Statistics Denmark and the Danish Health Data Authority. The data are available from the authorities, but resctrictions apply.

contact: Theresa Wimberley tw.ncrr@au.dk

project link: 
