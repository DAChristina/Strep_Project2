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
library(lubridate)

wd = "C:/Users/dac23/Downloads"
wd = "/home/ron/Downloads"
setwd(wd)

# Set dummy data
# SOURCE: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6128537/
dat <- read_excel("pone.0203205.s001.xlsx") %>% 
  glimpse()

dat <- dat %>% 
  rename(Week = `Epi Week`,
         Age_group = `Age group`) %>% 
  glimpse()

dates <- seq.Date(as.Date(paste0(2015, "-12-01")),
                  by = "week",
                  length.out = 70)

new_dat <- data.frame(
  Year = year(dates),
  Week = week(dates),
  Day = dates)

glimpse(new_dat)

new_datWeek <- left_join(dat, new_dat, by = c("Year", "Week"))

# Filter to Strep PN only (new_datWeek is a mixture of HI, NM, other Streps)
# See Etiology = SPN coz' others are mixture (HI, SPN, HI, etc), total should be 153
dat_SPN <- new_datWeek %>% 
  filter(Etiology == "SPN") %>% 
  glimpse()

# Specify to ST1 (total should be 85 cases)
dat_ST1 <- new_datWeek %>% 
  filter(PCR_SpSerotype == "1") %>% 
  glimpse()


# 1. Notes about data: #########################################################
# NO daily incidences, but weekly incidences available (constant aggregation)
# SOURCE: https://mrc-ide.github.io/EpiEstim/articles/EpiEstim_aggregated_data.html

# I have changed the data as originally as possible (1 row means 1 case)
# Crucial info:
# $ Day (also contains Year & Week)
# $ Age_group
# $ Pneumovac
# $ PCR_SpSerotype

# EpiEstim data required for estimate_R (based on > data("Flu2009")):
# incidence (n per-day)
# SI distribution (could be in separate df, data not available)

# Incidence viz per week (using library(incidence))
# SOURCE: https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html

# 1.1. (not) Directly plot the incidence
# SOURCE: https://www.repidemicsconsortium.org/incidence/
dat_ST1_i7 <- incidence(dat_ST1$Day, interval = 7)
# Unfortunately, 8 missing observations were removed (only year is known therefore Day is NA)

# By using library(incidence)
plot(dat_ST1_i7, # Or use dates = new_datWeek$Row_numb instead
     xlab = "Week number (1 Dec 2015 to 31 Mar 2017)",
     ylab = "Weekly Incidence of Serotype 1")

# 2. The EpiEstim trials #######################################################
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
output <- EpiEstim::estimate_R(incid = dat_ST1_i7$counts, # change df to new_datWeek or new_datWeek_ST1
                               dt = 7L,
                               dt_out = 7L,
                               iter = 10L, # 10L by default
                               tol = 1e-6, # 1e-6 by default
                               recon_opt = "match", # "naive" by default or "match"
                               method = "parametric_si",
                               config = config)

plot(output, legend = F,
     main = "The R(t) Estimation for Serotype 1 Incidence")


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


# Trial MCMC to generate SI from data cannot be executed because there are no si_data
# SOURCE: https://mrc-ide.github.io/EpiEstim/reference/Flu2009.html
if (FALSE) {
  ## Note the following examples use an MCMC routine
  ## to estimate the serial interval distribution from data,
  ## so they may take a few minutes to run
  
  ## estimate the reproduction number (method "si_from_data")
  output_MCMC <- EpiEstim::estimate_R(incid = new_datWeek$I_week,
                                      dt = 7L, # ???
                                      method="si_from_data",
                                      si_data = "???",
                                      config = make_config(
                                        list(mcmc_control = 
                                               make_mcmc_control(list(
                                                 burnin = 1000,
                                                 thin = 10,
                                                 seed = 1)),
                                             n1 = 1000, n2 = 50,
                                             si_parametric_distr = "G")))
  
  plot(output_MCMC, legend = F)
  ## the second plot produced shows, at each each day,
  ## the estimate of the reproduction number
  ## over the 7-day window finishing on that day.
}

# Example using SIR model but deterministic:
# https://www.kaggle.com/code/sunfinger/can-we-simply-predict-evolution-of-covid-19
