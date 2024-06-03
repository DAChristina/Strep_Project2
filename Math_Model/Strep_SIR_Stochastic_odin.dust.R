
# Running the SIR model with dust
pars <- list(A_ini = 6e7*(2e-6), # S_ini*(2e-6) = 120 people, 
             time_shift = 0.2,
             beta_0 = 0.06565,
             beta_1 = 0.07,
             wane = 0.002,
             log_delta = (-4.98), # will be fitted to logN(-10, 0.7)
             sigma_1 = (1/15.75),
             sigma_2 = (1)
) # Serotype 1 is categorised to have the lowest carriage duration

sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = 15L,
                         n_threads = 4L,
                         seed = 1L)

sir_model$state()

# update_state is required "every single time" to run & produce matrix output (don't know why)
sir_model$update_state(pars = pars,
                       time = 0) # make sure time is 0

# all_date <- incidence$day
# all_date <- data.frame(col = integer(4745))
n_times <- 4745 # 4745 or similar to the number of date range (of the provided data), or try 500 for trial
n_particles <- 15
x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

# Beta check
time <- seq(1, n_times, 1)
time_shift <- 70
beta <- pars$beta_0*(1+pars$beta_1*sin(2*pi*(time_shift+time)/365))
max(beta)
min(beta)

# R0 estimation (R0 changes due to seasonality)
R0 <- (beta/(pars$log_delta+pars$sigma_1)) +  ((pars$log_delta)*(beta)) / ((pars$log_delta + 192/(4064*4745))*(pars$sigma_2 + 192/(4064*4745))) # print R0
max(R0)
min(R0)
# plot(time, R0)
# pars$beta_1/(pars$delta) + (pars$qu*pars$delta)/(pars$delta*pars$sigma) # print R0


for (t in seq_len(n_times)) {
  x[ , , t] <- sir_model$run(t)
}
time <- x[1, 1, ] # because in the position of [1, 1, ] is time
x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
library(tidyverse)
glimpse(x)

## 1. Data Load ################################################################
library(readxl)
# New updated data with meningitis (25.04.2024)
dat <- read_excel("serotype1_UKHSA_imperial_date_age_region_MOLIS_withdeath_meningitis_clean.xlsx") #%>% 
# glimpse()

dat <- dat %>% 
  rename(Earliest.specimen.date = Earliestspecimendate,
         current.region.name = currentregionname)

dat_G <- dat %>% 
  mutate(AGEYR = ifelse(AGEYR >= 90, 90, as.numeric(AGEYR)), # For incidence calculation, data grouped for people aged 90+
         year = year(Earliest.specimen.date),
         month = month(Earliest.specimen.date),
         vacc = case_when(
           year < 2006 ~ "Pre-PCV7",
           year >= 2006 & year < 2011 ~ "PCV7",
           year >= 2011 ~ "PCV13",
           TRUE ~ NA_character_
         ),
         ageGroup = case_when(
           AGEYR < 2 ~ "<2",
           AGEYR >= 2 & AGEYR < 5 ~ "2-4",
           AGEYR >= 5 & AGEYR < 15 ~ "5-14",
           AGEYR >= 15 & AGEYR < 31 ~ "15-30", # Edit the Age-band into 15-30 & 31-44
           AGEYR >= 31 & AGEYR < 45 ~ "31-44", # Edit the Age-band into 15-30 & 31-44
           AGEYR >= 45 & AGEYR < 65 ~ "45-64",
           AGEYR >= 65 ~ "65+",
           is.na(AGEYR) ~ "Unknown" # 16 IDs have no AGEYR
           # TRUE ~ "Unknown" 
         ),
         current.region.name = ifelse(current.region.name == "EASTERN", "EAST", current.region.name), # Wrong perception of "EASTERN" that should be "EAST"
         current.region.name = case_when(
           current.region.name == "E MIDS" ~ "East Midlands",
           current.region.name == "EAST" ~ "East of England",
           current.region.name == "LONDON" ~ "London",
           current.region.name == "N EAST" ~ "North East",
           current.region.name == "N WEST" ~ "North West",
           current.region.name == "S EAST" ~ "South East",
           current.region.name == "S WEST" ~ "South West",
           current.region.name == "W MIDS" ~ "West Midlands",
           current.region.name == "YORK&HUM" ~ "Yorkshire and The Humber",
           TRUE ~ current.region.name
         ),
         ageLabel = ifelse(AGEYR >= 90, 90, as.numeric(AGEYR)), # For incidence calculation, data grouped for people aged 90+
  ) %>% 
  glimpse()

# Basic case count data without age structure or regions
# Create all hypothetical recorded disease date
dat_G$Earliest.specimen.date <- as.Date(dat_G$Earliest.specimen.date)
all_date <- data.frame(allDate = seq.Date(from = min(dat_G$Earliest.specimen.date),
                                          to = max(dat_G$Earliest.specimen.date), 
                                          by = 1))
all_date$day <- 1:nrow(all_date)
# Coz the incidence only requires 2 columns called "counts" and "Day" in NUMBERS
# The counts (but in 0 counts the date are not recorded)
Natm_ni <- dat_G %>% 
  group_by(Earliest.specimen.date) %>% 
  summarise(counts_Ser1 = n()) %>% 
  ungroup() #%>% 
# glimpse()

Natm_nmeningitis <- dat_G %>% 
  filter(MeningitisFlag == "Y") %>% 
  group_by(Earliest.specimen.date) %>% 
  summarise(counts_meningitis = n()) %>% 
  ungroup() #%>% 
# glimpse()

Natm_n30DDeath <- dat_G %>% 
  filter(`30daydeath` == "D") %>% 
  group_by(Earliest.specimen.date) %>% 
  summarise(counts_30DDeath = n()) %>% 
  ungroup() #%>% 
# glimpse()

# Create a new df based on counts per day for Serotype 1, meningitis, and 30 days death
Natm_n_i <- full_join(all_date, Natm_ni,
                      by = c("allDate" = "Earliest.specimen.date"))

Natm_n_im <- full_join(Natm_n_i, Natm_nmeningitis,
                       by = c("allDate" = "Earliest.specimen.date"))

Natm_n_imD <- full_join(Natm_n_im, Natm_n30DDeath,
                        by = c("allDate" = "Earliest.specimen.date")) %>% 
  replace(is.na(.), 0) #%>% # NA means no data of meningitis or 30 days death, changed them to 0
  #glimpse()

# Total population data by age, year for each region
# SOURCE: https://www.nomisweb.co.uk/
# pop <- read_excel("nomis_2024_04_15_124553_DCedit.xlsx") %>% 
# glimpse()
# I don't think I need total population for now,
# Examples on https://github.com/mrc-ide/mcstate/blob/master/inst/sir_incidence.csv
# Requires case count per aligned day only

# Viz per-day counts by base R plot
par(bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))

col_imD <- c(counts_Ser1 = "deepskyblue3",
             counts_meningitis = "green",
             counts_30DDeath = "maroon")


plot(Natm_n_imD$allDate, Natm_n_imD$counts_Ser1, type = "b",
     xlab = "Date", ylab = "Counts",
     ylim = c(0, max(Natm_n_imD$counts_Ser1)+2),
     col = col_imD[1], pch = 20)

lines(Natm_n_imD$allDate, Natm_n_imD$counts_meningitis,
      type = "b", col = col_imD[2], pch = 20)
lines(Natm_n_imD$allDate, Natm_n_imD$counts_30DDeath,
      type = "b", col = col_imD[3], pch = 20)
legend("topleft", names(col_imD), fill = col_imD, bty = "n")

## 2. Data Fitting #############################################################
# The anatomy of an mcstate particle filter, as noted above, consists of three main components: \n 
# 1. A set of observations to fit the model to, generated using mcstate::particle_filter_data(). \n 
# 2. A model to fit, which must be a dust generator, either dust::dust() or odin.dust::odin_dust(). \n 
# 3. A comparison function, which is an R function which calculates the likelihood of the state given the data at one time point.


# There is a calibration function in mcstate to fit our model to data.
# https://mrc-ide.github.io/mcstate/articles/sir_models.html

incidence <- Natm_n_imD %>% 
  mutate(day_calibrated = day+(4745)) %>% # day_calibrated to simulate disease epidemic happens 13 years after the outbreak (4745*2)
  select(day, day_calibrated, counts_Ser1) %>% 
  rename(cases = counts_Ser1) # That annoying name

par(mar = c(5.1, 5.1, 0.5, 0.5), mgp = c(3.5, 1, 0), las = 1)
cols <- c(S = "#8c8cd9", A = "darkred", D = "#cc0099", R = "#999966", n_AD_daily = "orange", n_AD_cumul = "green")
# matplot(time, t(x[1, , ]), type = "l",
#         xlab = "Time", ylab = "Number of individuals",
#         col = cols[["S"]], lty = 1, ylim = range(x))
matplot(time, t(x[3, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["D"]], lty = 1)#, ylim = max(x[2,,]))

matlines(incidence$day, incidence$cases, type = "l", col = "steelblue")

# matlines(time, t(x[2, , ]), col = cols[["A"]], lty = 1)
# matlines(time, t(x[3, , ]), col = cols[["D"]], lty = 1)
# matlines(time, t(x[4, , ]), col = cols[["R"]], lty = 1)
# matlines(time, t(x[5, , ]), col = cols[["n_AD_daily"]], lty = 1)
# matlines(time, t(x[6, , ]), col = cols[["n_AD_cumul"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
max(x[5,,]) # Check max n_AD_daily
max(x[3,,]) # Check max D

# write.csv(x, file="Output_sir_result.csv", row.names =T)
