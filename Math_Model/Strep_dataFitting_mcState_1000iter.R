
## 1. Data Preparation #########################################################
if (!require(mcstate, quietly=T)){
  install.packages("mcstate",
                   repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
  
  library(mcstate)
}

if (!require(odin.dust, quietly=T)){
  install.packages("odin.dust",
                   repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
  
  library(odin.dust)
}


library(tidyverse)
library(readxl)
library(coda)

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
  replace(is.na(.), 0) %>% # NA means no data of meningitis or 30 days death, changed them to 0
  glimpse()

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
  select(day, counts_Ser1) %>% 
  rename(cases = counts_Ser1) # That annoying name
# To make my life easier I compile the Serotype 1 cases into a new object called sir_data
# data is fed as an input to mcstate::particle_filter_data
dt <- 1 # rate must be an integer; 0.25 to make it 4 days, I make it 1
sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "day",
                                          rate = 1 / dt)
rmarkdown::paged_table(sir_data)
# And ofc to make sure this is Strep's data:
plot(incidence$day, incidence$cases,
     type = "l", xlab = "Day", ylab = "New cases")


## 2a. Model Load ##############################################################
# The model below is stochastic, closed system SADR model that I have created before
# I update the code, filled the mprameters with numbers;
# e.g.dt <- user(1) because if dt <- user() generates error during MCMC run
gen_sir <- odin.dust::odin_dust("sir_stochastic.R")

# This is part of sir odin model:
pars <- list(dt = 1,
             S_ini = 6e7, # England's pop size is roughly 67,000,000
             A_ini = 100,
             D_ini = 0,
             time_shift = 71.88781655,
             beta_0 = 0.0645,
             beta_1 = 0.07,
             log_delta = (-4),
             sigma_1 = (1/15.75), # FIXED carriage duration of diseased = 15.75 days (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
             sigma_2 = (1) # FIXED estimated as acute phase
) # Serotype 1 is categorised to have the lowest carriage duration

gen_sir$new(pars = pars,
            time = 1,
            n_particles = 10L,
            n_threads = 4L,
            seed = 1L)$info()

# https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model-over-time
mod <- gen_sir$new(pars, 0, 20)
y <- mod$simulate(c(0, sir_data$time_end))
i <- mod$info()$index[["time"]]
j <- mod$info()$index[["n_AD_daily"]]
matplot(y[i, 1, ], t(y[j, , ]), type = "l", col = "#00000055", lty = 1, las = 1,
        xlab = "Day", ylab = "Cases")
points(cases ~ day, incidence, col = "red", pch = 19)

index <- function(info) {
  list(run = c(incidence = info$index$n_AD_daily),
       state = c(t = info$index$time,
                 D = info$index$D,
                 n_AD_daily = info$index$n_AD_daily))
}
index(mod$info())


## 2b. The Comparison Function #################################################
# Further details: https://mrc-ide.github.io/mcstate/articles/sir_models.html
# gen_sir$new(pars = list(), time = 0, n_particles = 1L)$info()

case_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6
  
  incidence_modelled <- state[6, , drop = TRUE] # (incidence based on model "n_AD_daily")
  incidence_observed <- observed$cases # daily new cases
  lamb <- incidence_modelled +
    rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lamb, log = TRUE)
}

# If use SIR example the calculation below is not required:
# incidence_compare <- function(state, prev_state, observed, pars = NULL) {
#   exp_noise <- 1e6
# 
#   lambda <- state[7, , drop = TRUE] +
#     rexp(n = length(incidence_modelled), rate = exp_noise) # lambda = n_SI_cumul + rexp(incid, rate)
#   dpois(x = observed$cases, lambda = lambda, log = TRUE)
# }

# Plot these along with the data
# true_history seems like a simulated data consisting of SIR plus cases
# see: https://mrc-ide.github.io/mcstate/articles/restart.html

# Imma create the file in odin based on beta and sigma, a basic SIR model, closed system,
# with number of rows equivalent to the current output, day = nrow(all_date)
# transpose is not required because this is the output of sir.odin
sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = 1L,
                         n_threads = 4L,
                         seed = 1L)

# update_state is required "every single time" to run & produce matrix output (don't know why)
sir_model$update_state(pars = pars,
                       time = 0) # make sure time is 0

all_date <- incidence$day
n_times <- length(all_date) # 4745 or similar to the number of date range (of the provided data), or try 500 for trial
n_particles <- 1
sir_output <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

for (t in seq_len(n_times)) {
  sir_output[ , , t] <- sir_model$run(t)
}
time <- sir_output[1, 1, ] # because in the position of [1, 1, ] is time
sir_output <- sir_output[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
glimpse(sir_output)


# weird things happen here. true_history example:
# true_history <- readRDS("sir_true_history.rds")
# seems like they combine the SIR model output with n_SI_daily from the model,
# make them as a wide-type arrays with 4 rows (S,I,R,n_SI_daily)
# create a new matrix with initial state (t1 = 0)
begin_t0_value <- c(6e7, 0, 0, 0, 0, 0)  # Receall S_ini in Pars
# begin_t0 <- array(begin_t0_value, dim = c(5, 1, 1))
begin_t0 <- matrix(nrow = 6, ncol = 1, data = begin_t0_value)

# The sir_output with t1 = 1
chosen_output <- sir_output[1:6, ]
bindedMtx <- as.matrix(chosen_output)

# Combine
true_history <- cbind(begin_t0, bindedMtx) # create a 2D matrix and then change them into 3D:
true_history <- array(bindedMtx, dim = c(6, 1, length(all_date)+1))
glimpse(true_history)
# ori_history <- array(bindedMtx, dim = c(4, 1, length(all_date))) # 1-4 means S, I, R, n_SI_daily
# true_history <- array(c(begin_t0, ori_history), dim = c(dim(begin_t0)[1], dim(ori_history)[2], length(all_date)+1)) # 3-d matrices bind




# recall true_history
plot_particle_filter <- function(history, true_history, times, obs_end = NULL) {
  if (is.null(obs_end)) {
    obs_end <- max(times)
  }
  
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  cols <- c(S = "#8c8cd9", A = "darkred", D = "#cc0099", R = "#999966", n_AD_daily = "orange", n_AD_cumul = "green")
  matplot(times, t(history[4, , -1]), type = "l", # I change history[2, , -1] becaue history[1, , -1] is time
          xlab = "Time", ylab = "Number of individuals",
          col = cols[["D"]], lty = 1) #, ylim = range(history))
  # matlines(times, t(history[3, , -1]), col = cols[["A"]], lty = 1)
  # matlines(times, t(history[4, , -1]), col = cols[["D"]], lty = 1)
  # matlines(times, t(history[5, , -1]), col = cols[["R"]], lty = 1)
  matpoints(times[1:obs_end], t(true_history[1:3, , -1]), pch = 1,
            col = "green")
  legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
}

# Inferring Parameters
n_particles <- 500 # I increase the particles from 100 to 500
filter <- mcstate::particle_filter$new(data = sir_data,
                                       model = gen_sir, # Changing gen_sir into sir_model (with updated parameters that I can control) produce error:
                                       # Error in initialize(...) : 'model' must be a dust_generator
                                       n_particles = n_particles,
                                       compare = case_compare,
                                       seed = 1L)

# recall pars but dt = dt, and dt <- 1
dt <- 1
filter$run(save_history = TRUE, pars = list(dt = 1,
                                            S_ini = 6e7, # England's pop size is roughly 67,000,000
                                            A_ini = 100,
                                            D_ini = 0,
                                            time_shift = 71.88781655,
                                            beta_0 = 0.0645,
                                            beta_1 = 0.07,
                                            log_delta = (-4),
                                            sigma_1 = (1/15.75), # carriage duration of diseased = 15.75 days (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
                                            sigma_2 = (1) # estimated as acute phase
) # Serotype 1 is categorised to have the lowest carriage duration
)

plot_particle_filter(filter$history(), true_history, incidence$day)


## 3. MCMC Run #################################################################
# Using MCMC to Infer Parameters (Metropolis-Hastings):
# https://mrc-ide.github.io/mcstate/reference/pmcmc_parameters.html
# first time I run the code & it requires ~10' for 500 steps, 200 burnin

# Invasiveness estimate based on https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009389
# S1 Dataset. Invasiveness estimates for serotypes in children.
# Serotype	invasiveness	invasiveness_lower	invasiveness_upper
# 1	0.017369780939027	0.009888269105106	0.031491598702714

# S2 Dataset. Invasiveness estimates for serotypes in adults.
# Serotype	invasiveness	invasiveness_lower	invasiveness_upper
# 1	0.019837434032008	0.007774465457634	0.057617171798831
# 2 Dataset. Invasiveness estimates for serotypes in adults.
# Serotype	
# For this trial I use the mean number of invasiveness estimate:
time_shift <- mcstate::pmcmc_parameter("time_shift", 71.88781655, min = 0, prior = function(s)
  dunif(s, min = 0, max = 365/2, log = TRUE)) # ~Uniform[0,365/2]
beta_0 <- mcstate::pmcmc_parameter("beta_0", 0.0645, min = 0, prior = function(q)
  dgamma(q, shape = 1, scale = 0.1, log = TRUE)) # draws from gamma distribution dgamma(1, 0.2) --> exp dist
beta_1 <- mcstate::pmcmc_parameter("beta_1", 0.07, min = 0, prior = function(r)
  dgamma(r, shape = 1, scale = 0.1, log = TRUE)) # draws from gamma distribution dgamma(1, 0.2) --> exp dist
# For dGamma, I change:
# shape into prior mean^2/variance, given prior mean = init, variance = 0.1 (larger, instead of 0.01)
# scale into prior mean/variance, given prior mean = init, variance = 0.1 (larger, instead of 0.01)
# beta_0 <- mcstate::pmcmc_parameter("beta_0", 0.0645, min = 0, prior = function(q)
#   dgamma(q, shape = (((pars$beta_0)^2)/0.1), scale = (((pars$beta_0))/0.1), log = TRUE)) # draws from gamma distribution dgamma(1, 0.2) --> exp dist
# beta_1 <- mcstate::pmcmc_parameter("beta_1", 0.07, min = 0, prior = function(r)
#   dgamma(r, shape = (((pars$beta_1)^2)/0.1), scale = (((pars$beta_1))/0.1), log = TRUE)) # draws from gamma distribution dgamma(1, 0.2) --> exp dist

log_delta <- mcstate::pmcmc_parameter("log_delta", (-4), prior = function(p)
  dunif(p, min = -5, max = 0.7, log = TRUE)) # logN distribution for children & adults (Lochen et al., 2022)

proposal_matrix <- diag(0.1, 4) # assumption no co-variance occur
# proposal_matrix <- as.matrix(0.1)
mcmc_pars <- mcstate::pmcmc_parameters$new(list(time_shift = time_shift,
                                                beta_0 = beta_0,
                                                beta_1 = beta_1,
                                                log_delta = log_delta),
                                           proposal_matrix)


n_steps <- 500
n_burnin <- n_steps/2

control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE,
  save_trajectories = TRUE,
  # rerun_every = 7,
  progress = TRUE)
pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)
plot_particle_filter(pmcmc_run$trajectories$state, true_history, incidence$day)

processed_chains <- mcstate::pmcmc_thin(pmcmc_run, burnin = n_burnin, thin = 2)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
parameter_mean_hpd

mcmc1 <- coda::as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))
summary(mcmc1)
plot(mcmc1)


# Diagnostics TRIAL ############################################################
# coz usually the first run is not efficient (in terms of effective sample size per iteration)
print("pMCMC Effective Size & Rejection Rate")
coda::effectiveSize(mcmc1)
1 - coda::rejectionRate(mcmc1)

# Autocorrelation plots
print("Autocorrelation?")
AuCorr_time_shift <- coda::acfplot(mcmc1[, "time_shift"], main = "Autocorrelation time shift")
AuCorr_time_shift
AuCorr_beta_0 <- coda::acfplot(mcmc1[, "beta_0"], main = "Autocorrelation beta0")
AuCorr_beta_0
AuCorr_beta_1 <- coda::acfplot(mcmc1[, "beta_1"], main = "Autocorrelation beta1")
AuCorr_beta_1
AuCorr_log_delta <- coda::acfplot(mcmc1[, "log_delta"], main = "Autocorrelation log(delta)")
AuCorr_log_delta
## 3a. Tuning the pMCMC part 1 #################################################
# Use the covariance of the state as the proposal matrix:
proposal_matrix <- cov(pmcmc_run$pars)
mcmc_pars <- mcstate::pmcmc_parameters$new(
  list(time_shift = time_shift,
       beta_0 = beta_0,
       beta_1 = beta_1,
       log_delta = log_delta),
  proposal_matrix)
proposal_matrix

control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE,
  save_trajectories = TRUE,
  # rerun_every = 7,
  progress = TRUE,
  n_chains = 4)
pmcmc_tuned_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)

mcmc2 <- coda::as.mcmc(cbind(
  pmcmc_tuned_run$probabilities, pmcmc_tuned_run$pars))

summary(mcmc2)
plot(mcmc2)

print("Tuning Result Effective Size & Rejection Rate")
coda::effectiveSize(mcmc2)
1 - coda::rejectionRate(mcmc2)

# Gelman-Rubin diagnostic
# https://cran.r-project.org/web/packages/coda/coda.pdf
# Extract the number of chains
n_chains <- control$n_chains

# Determine the length of each chain
n_samples <- nrow(pmcmc_tuned_run$pars) / n_chains

# Split the parameter samples and probabilities by chains
chains <- lapply(1:n_chains, function(i) {
  start <- (i - 1) * n_samples + 1
  end <- i * n_samples
  list(
    pars = pmcmc_tuned_run$pars[start:end, ],
    probabilities = pmcmc_tuned_run$probabilities[start:end, ]
  )
})

# Convert chains to mcmc objects
mcmc_chains <- lapply(chains, function(chain) {
  as.mcmc(cbind(chain$probabilities, chain$pars))
})

# Combine chains into a list
mcmc_chains_list <- do.call(mcmc.list, mcmc_chains)

print("Gelman-Rubin diagnostic")
coda::gelman.plot(mcmc_chains_list,
                  bin.width = 10,
                  max.bins = 50,
                  confidence = 0.95,
                  transform = FALSE,
                  autoburnin=TRUE,
                  auto.layout = TRUE)
# ask, col, lty, xlab, ylab, type, ...)

coda::gelman.diag(mcmc_chains_list,
                  confidence = 0.95,
                  transform=FALSE,
                  autoburnin=TRUE,
                  multivariate=F) # Change multivariate = F instead of T

## 4. Running predictions
