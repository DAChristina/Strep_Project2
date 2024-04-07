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

dates <- seq.Date(as.Date(paste0(2015, "-01-01")),
                  as.Date(paste0(2017, "-12-31")),
                  by = "week")

new_dat <- data.frame(
  Year = year(dates),
  Week = week(dates)
) %>% 
  glimpse()

new_datWeek <- left_join(new_dat, dat_week, by = c("Year", "Week"))
new_datWeek <- new_datWeek %>% 
  mutate(I_Week = if_else(is.na(I_week), 0, I_week)) %>% 
  filter(I_week != 0) %>% 
  glimpse()


# Incidence viz per week (using library(incidence))
# SOURCE: https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
plot(as.incidence(new_datWeek$I_week, dates = new_datWeek$Week))

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

# Example using SIR model but deterministic:
# https://www.kaggle.com/code/sunfinger/can-we-simply-predict-evolution-of-covid-19
