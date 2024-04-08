# EpiEstim & projections are part of RECON:
# https://www.repidemicsconsortium.org/projects/

if(!require(EpiEstim, quietly=T)){
  # Install the package from repos to get an updated version of EpiEstim!
  # Dependencies required in deb: sudo apt-get install libsodium-dev
  install.packages('EpiEstim', repos = c('https://mrc-ide.r-universe.dev', 'https://cloud.r-project.org'))
  library(EpiEstim)
}

library(tidyverse)
library(readxl)
library(incidence)

wd = "C:/Users/dac23/Downloads"
wd = "/home/ron/Downloads"
setwd(wd)

# Set dummy data
# SOURCE: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6128537/
dat <- read_excel("pone.0203205.s001.xlsx") %>% 
  glimpse()

sort(unique(dat$Latex))
sort(unique(dat$Culture))
sort(unique(dat$PCR))
sort(unique(dat$Etiology))
# See Etiology = SPN coz' others are mixture (HI, SPN, HI, etc), total should be 153
sort(unique(dat$PCR_SpSerotype))

dat <- dat %>% 
  filter(Etiology == "SPN") %>% 
  rename(Week = `Epi Week`,
         Age_group = `Age group`) %>% 
  glimpse()

# 1. Notes about data: #########################################################
# NO daily incidences, but weekly incidences available (constant aggregation)
# SOURCE: https://mrc-ide.github.io/EpiEstim/articles/EpiEstim_aggregated_data.html

# Crucial info:
# $ Year
# $ Week
# $ Age_group
# $ Pneumovac
# $ PCR_SpSerotype

# EpiEstim data required for estimate_R (based on > data("Flu2009")):
# incidence (n per-day)
# SI distribution (could be in separate df)

# Incidence count per week
dat_week <- dat %>% 
  group_by(Year, Week) %>% 
  summarise(I_week = n()) %>% 
  ungroup() %>% 
  arrange(Year, Week) %>% 
  glimpse()

# Incidence count per week for Serotype 1
dat_week_ST1 <- dat %>% 
  filter(PCR_SpSerotype == "1") %>% 
  group_by(Year, Week) %>% 
  summarise(I_week = n()) %>% 
  ungroup() %>% 
  # filter(!is.na(Week)) %>% # filtered because week unknown
  # mutate(Week = if_else(is.na(Week), 0, Week)) %>% 
  glimpse()


# Let's arrange the weeks
library(lubridate)

dates <- seq.Date(as.Date(paste0(2015, "-12-01")),
                  by = "week",
                  length.out = 70)

new_dat <- data.frame(
  Year = year(dates),
  Week = week(dates),
  Day = dates)

glimpse(new_dat)

new_datWeek <- left_join(new_dat, dat_week, by = c("Year", "Week"))
new_datWeek <- new_datWeek %>% 
  mutate(I_week = if_else(is.na(I_week), 0, as.numeric(I_week)),
         Row_numb = row_number()) %>% 
  # filter(I_week != 0) %>% 
  glimpse()


# Incidence viz per week (using library(incidence))
# SOURCE: https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html

# 1.1. Directly plot the incidence
plot(as.incidence(new_datWeek$I_week, dates = new_datWeek$Day), # Or use dates = new_datWeek$Row_numb instead
     xlab = "Week number (1 Dec 2015 to 31 Mar 2017)",
     ylab = "Weekly incidence")

# 1.2. Or save the data as incidence (FAILED)
SPN_incid <- as.incidence(new_datWeek$I_week, dates = new_datWeek$Day)

plot(SPN_incid) %>%
  add_projections(proj, c(0.025, 0.5, 0.975))

# 2. The trials ################################################################
# SOURCE: https://mrc-ide.github.io/EpiEstim/articles/EpiEstim_aggregated_data.html

# Serial interval (SI):
# Average time between the symptoms of a primary case to a secondary case.
# Carriage duration:
# A period when an individual carries and spreads the infectious agents

# As far as I know, carriage duration is not always equivalent to serial interval.
# For this trial, assumed that carriage duration is equivalent to SI

# Literature to specify a parametric SI (Chaguza et al., 2021):
# mean_SI = unknown, assumed as 15.75 days (95% CI 7.88-31.49) (Serotype 1)
# SD_SI = unknown, assumed as sqrt(8)*(31.49-7.88)/(2*1.96)
# Dist = unknown

mean_SI = 15.75
SD_SI = sqrt(8)*(31.49-7.88)/(2*1.96) # 17.0355
config = make_config(list(mean_si = mean_SI,
                          std_si = SD_SI))

# EpiEstim but incidence should be in days (not weeks):
# output <- EpiEstim::estimate_R(incid = dat_week$I_week,
# method = "parametric_si",
# config = config)

# Previously failed, estimate_R can run for weekly data:
output <- EpiEstim::estimate_R(incid = new_datWeek$I_week, # change df to dat_week or dat_week_ST1
                               dt = 7L,
                               dt_out = 7L,
                               iter = 10L, # 10L by default
                               tol = 1e-6, # 1e-6 by default
                               recon_opt = "match", # "naive" by default or "match"
                               method = "parametric_si",
                               config = config)

plot(output, legend = F)


# Forecast trial (data should be daily incidence)
# SOURCE: https://mrc-ide.github.io/EpiEstim/articles/full_EpiEstim_vignette.html
if(!require(projections, quietly=T)){
  devtools::install_github("reconhub/projections")
  library(projections)
}

trunc_date <- max(new_datWeek$Day) - 1 # because data report already in weeks
trunc_linelist <- subset(new_datWeek, new_datWeek$Day < trunc_date) # counts = I_week = Incidence counts

R_si_parametric_recent <- estimate_R(incid = trunc_linelist$I_week, 
                                     method = "parametric_si",
                                     config = make_config(mean_si = mean_SI,
                                                          std_si = SD_SI,
                                                          t_start = length(trunc_linelist$I_week) - 1,
                                                          t_end = length(trunc_linelist$I_week)))



# trunc_linelist convert from weekly to daily:
trunc_linelist_daily <- trunc_linelist %>%
  complete(Day, nesting(Week)) %>%
  select(-Week) %>%
  mutate(I_day = I_week / 7,
         I_day = if_else(is.na(I_day), 0, as.numeric(I_day))) %>% 
  glimpse()

evd_incid_trunc <- as.incidence(trunc_linelist_daily$I_day, dates = trunc_linelist_daily$Day)

# Ofc FAILED because project() requires daily incidence data
# Eventhouth the weekly incidence data was converted to daily incidence,
# It doesn't work because of the data type in Incidence is decimal (not integer)
proj <- project(evd_incid_trunc, # truncated incidence object
                R = R_si_parametric_recent$R$`Median(R)`, # R estimate
                si = R_si_parametric_recent$si_distr[-1], # SI (starting on day 1)
                n_sim = 1000, # simulate 1000 trajectories
                n_days = 52, # 52 weeks = a year
                R_fix_within = TRUE) # keep the same value of R every day

plot(as.incidence(new_datWeek$I_week, dates = new_datWeek$Day), # Or use dates = new_datWeek$Row_numb instead
     xlab = "Week number (1 Dec 2015 to 31 Mar 2017)",
     ylab = "Weekly incidence") %>%
  add_projections(proj, c(0.025, 0.5, 0.975))



# Example using SIR model but deterministic:
# https://www.kaggle.com/code/sunfinger/can-we-simply-predict-evolution-of-covid-19
