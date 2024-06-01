
## 1. Data Preparation #########################################################
library(mcstate)
library(odin.dust)

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

hist(incidence$cases) # huge zero daily cases occur

# To make my life easier I compile the Serotype 1 cases into a new object called sir_data
# data is fed as an input to mcstate::particle_filter_data
dt <- 1 # rate must be an integer; 0.25 to make it 4 days, I make it 1
sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "day",
                                          rate = 1 / dt)

# Annotate the data so that it is suitable for the particle filter to use
dat <- mcstate::particle_filter_data(incidence, "day", 4, 0)

rmarkdown::paged_table(sir_data)
# And ofc to make sure this is Strep's data:
plot(incidence$day, incidence$cases, col = "deepskyblue3",
     type = "l", xlab = "Day", ylab = "New cases")


## 2a. Model Load ##############################################################
# The model below is stochastic, closed system SADR model that I have created before
# I update the code, filled the mprameters with numbers;
# e.g.dt <- user(1) because if dt <- user() generates error during MCMC run
gen_sir <- odin.dust::odin_dust("sir_stochastic.R")

# This is part of sir odin model:
pars <- list(time_shift = 72,
             beta_0 = 0.06565,
             beta_1 = 0.07,
             wane = 0.002,
             log_delta = (-4.98), # will be fitted to logN(-7, 0.7)
             sigma_2 = 1
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
matplot(y[i, 1, ], t(y[j, , ]), type = "l", col = "maroon", lty = 1, las = 1,
        xlab = "Day", ylab = "Cases", ylim = c(0, 10))
points(cases ~ day, incidence, col = "deepskyblue3", pch = 19)

index <- function(info) {
  list(run = c(incidence = info$index$n_AD_daily),
       state = c(t = info$index$time,
                 D = info$index$D,
                 n_AD_daily = info$index$n_AD_daily))
}
index(mod$info())


## 2b. The Comparison Function #################################################
# Further details: https://mrc-ide.github.io/mcstate/articles/sir_models.html
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model
# SEE the location of n_AD_daily from:
# gen_sir$new(pars = list(), time = 0, n_particles = 1L)$info()
# (ABOVE)

case_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e4
  
  incidence_modelled <- state[6, , drop = TRUE] # (incidence based on model's "n_AD_daily" from gen_sir)
  incidence_observed <- observed$cases # daily new cases
  lamb <- incidence_modelled +
    rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lamb, log = TRUE)
}

# That transform function
# https://github.com/mrc-ide/mcstate/blob/da9f79e4b5dd421fd2e26b8b3d55c78735a29c27/tests/testthat/test-if2.R#L40
# https://github.com/mrc-ide/mcstate/issues/184
parameter_transform <- function(pars) {
  time_shift <- pars[["time_shift"]]
  beta_0 <- pars[["beta_0"]]
  beta_1 <- pars[["beta_1"]]
  wane <- pars[["wane"]]
  log_delta <- pars[["log_delta"]]
  sigma_2 <- pars[["sigma_2"]]
  
  list(time_shift = time_shift,
       beta_0 = beta_0,
       beta_1 = beta_1,
       wane = wane,
       log_delta = log_delta,
       sigma_2 = sigma_2)
}

transform <- function(pars) {
  parameter_transform(pars)
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
# true_history seems like a simulated data consisting of SIR plus daily cases, IN 2d MATRIX FORMAT
# n_particles = 1
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
n_particles <- 1 # n_particles refers to n_particles in sir_model (=1)
sir_output <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

for (t in seq_len(n_times)) { # seq(0:seq_len(n_times), 1) failed to run
  sir_output[ , , t] <- sir_model$run(t)
}
time <- sir_output[1, 1, ] # because in the position of [1, 1, ] is time
# sir_output <- sir_output[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
glimpse(sir_output)


# weird things happen here. true_history example:
# true_history <- readRDS("sir_true_history.rds")
# seems like they combine the SIR model output with n_SI_daily from the model,
# make them as a wide-type arrays with 4 rows (S,I,R,n_SI_daily)
# create a new matrix with initial state (t1 = 0)
begin_t0_value <- c(6e7, 0, 0, 0, 0, 0, 0)  # Recall S_ini in Pars
# begin_t0 <- matrix(nrow = 7, ncol = 1, data = begin_t0_value)
begin_t0 <- array(begin_t0_value, dim=c(7, 10, 1)) # check dim in glimpse(sir_output)

# sir_output with t1 = 1 as a matrix
chosen_output <- sir_output[1:7, , ]
bindedMtx <- as.matrix(chosen_output)
glimpse(bindedMtx)

# Combine
true_history <- cbind(begin_t0, bindedMtx) # create a 2D matrix and then change them into 3D:
true_history <- array(bindedMtx, dim = c(7, 1, length(all_date)+1))
glimpse(true_history)


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
  matpoints(times[1:obs_end], t(true_history[2:4, , -1]), pch = 1,
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

filter$run(save_history = TRUE, pars = list(time_shift = 72,
                                            beta_0 = 0.06565,
                                            beta_1 = 0.07,
                                            wane = 0.002,
                                            log_delta = (-4.98) # will be fitted to logN(-5, 0.7)
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
# For this trial I use the mean number of invasiveness estimate: ((0.017369780939027+0.019837434032008)/2)/365 = 5.096879e-05 (Africa)
# Multiplying hypothetical delta by UK_calibration: 0.8066608*5.096879e-05 = 4.111452e-05
prepare_parameters <- function(initial_pars, priors, proposal, transform) {
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(
    list(mcstate::pmcmc_parameter("time_shift", 72.5, min = 1, max = 365,
                                  prior = function(s) dunif(s, min = 1, max = 365, log = TRUE)), # ~Uniform[72,73]
         mcstate::pmcmc_parameter("beta_0", 0.06565, min = 0, max = 0.8,
                                  prior = function(s) dgamma(s, shape = 1, scale = 0.1, log = TRUE)), # draws from gamma distribution dgamma(1, 0.2) --> exp dist)
         mcstate::pmcmc_parameter("beta_1", 0.07, min = 0, max = 0.8,
                                  prior = function(s) dgamma(s, shape = 1, scale = 0.1, log = TRUE)), # draws from gamma distribution dgamma(1, 0.2) --> exp dist
         mcstate::pmcmc_parameter("wane", 0.002, min = 0, max = 0.8,
                                  prior = function(s) dgamma(s, shape = 1, scale = 0.1, log = TRUE)), # draws from gamma distribution dgamma(1, 0.2) --> exp dist
         mcstate::pmcmc_parameter("log_delta", (-4.7), min = (-5), max = 0.7,
                                  prior = function(s) dunif(s, min = (-7), max = 0.7, log = TRUE)), # logN distribution for children & adults (Lochen et al., 2022)
         mcstate::pmcmc_parameter("sigma_2", 1, min = 0, max = 10,
                                  prior = function(s) dgamma(s, shape = 1, scale = 1, log = TRUE)) # shape = 1 , scale = 1 to capture 5 days/more (or dgamma(2.5, 0.5))?
    ),
    proposal = proposal,
    transform = transform)
  
}

prepare_priors <- function(pars){
  priors <- list()
  
  priors$time_shift <- function(s) {
    dunif(s, min = 1, max = 365, log = TRUE)
  }
  priors$beta_0 <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$beta_1 <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$wane <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$log_delta <- function(s) {
    dunif(s, min = (-7), max = 0.7, log = TRUE)
  }
  priors$sigma_2 <- function(s) {
    dgamma(s, shape = 1, scale = 1, log = TRUE)
  }
}

# Recall transform function
priors <- prepare_priors(pars)
# proposal_matrix <- diag(10, 6) # assumption no co-variance occur, variance in 10 days
proposal_matrix <- matrix(30, nrow = 6, ncol = 6)
rownames(proposal_matrix) <- c("time_shift", "beta_0", "beta_1", "wane", "log_delta", "sigma_2")
colnames(proposal_matrix) <- c("time_shift", "beta_0", "beta_1", "wane", "log_delta", "sigma_2")

mcmc_pars <- prepare_parameters(initial_pars = pars, priors = priors, proposal = proposal_matrix, transform = transform)

n_steps <- 1000
n_burnin <- n_steps/2

control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE,
  save_trajectories = TRUE,
  # rerun_every = 50,
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
print("pMCMC Effective Size & Acceptance Rate")
coda::effectiveSize(mcmc1)
1 - coda::rejectionRate(mcmc1)

# Autocorrelation plots
# print("Autocorrelation of mcmc1?")

# coda::acfplot(mcmc1[, "time_shift"], main = "Autocorrelation time shift")
coda::acfplot(mcmc1[, "beta_0"], main = "Autocorrelation beta0")
coda::acfplot(mcmc1[, "beta_1"], main = "Autocorrelation beta1")
coda::acfplot(mcmc1[, "wane"], main = "Autocorrelation wane")
coda::acfplot(mcmc1[, "log_delta"], main = "Autocorrelation log(delta)")
coda::acfplot(mcmc1[, "sigma_2"], main = "Autocorrelation sigma2")


## 3a. Tuning the pMCMC part 1 #################################################
# Use the covariance of the state as the proposal matrix:
proposal_matrix <- cov(pmcmc_run$pars)
mcmc_pars <- prepare_parameters(initial_pars = pars, priors = priors, proposal = proposal_matrix, transform = transform)

control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE,
  save_trajectories = TRUE,
  # rerun_every = 50,
  progress = TRUE,
  n_chains = 4)
pmcmc_tuned_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)

mcmc2 <- coda::as.mcmc(cbind(
  pmcmc_tuned_run$probabilities, pmcmc_tuned_run$pars))

summary(mcmc2)
plot(mcmc2)

print("Tuning Result Effective Size & Acceptance Rate of mcmc2")
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

print("Covariance matrix of mcmc2")
cov(as.matrix(mcmc_chains_list))

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


# Additional diagnostics #######################################################
# Source: https://mc-stan.org/bayesplot/reference/MCMC-scatterplots.html
library(bayesplot)
library(gridExtra)
# MCMC Pairs
bayesplot::mcmc_pairs(mcmc_chains_list,
                      pars = c("beta_0", "beta_1", "wane", "log_delta", "sigma_2"),
                      off_diag_args = list(size = 0.75))

bayesplot::mcmc_pairs(mcmc_chains_list,
                      off_diag_args = list(size = 0.75))

pairs(data.frame(mcmc2),
      main = "Another Pairs Visualisation")


# MCMC Scatter Pairs (I take only 3 parms; beta_0, beta_1, and log_delta)
par(mfrow = c(2,2))
color_scheme_set("orange")
beta_0_vs_beta_1 <- bayesplot::mcmc_scatter(
  mcmc_chains_list, 
  pars = c("beta_0", "beta_1"),
  size = 1.5,
  alpha = 0.25
)
beta_0_vs_beta_1 + stat_density_2d(color = "black", size = .5)

beta_0_vs_log_delta <- bayesplot::mcmc_scatter(
  mcmc_chains_list, 
  pars = c("beta_0", "log_delta"),
  size = 1.5,
  alpha = 0.25
)
beta_0_vs_log_delta + stat_density_2d(color = "black", size = .5)

beta_1_vs_log_delta <- bayesplot::mcmc_scatter(
  mcmc_chains_list, 
  pars = c("beta_1", "log_delta"),
  size = 1.5,
  alpha = 0.25
)
beta_1_vs_log_delta + stat_density_2d(color = "black", size = .5)

par(mfrow = c(1,1))


# MCMC Divergence (re-sampling by using pmcmc_sample)
# https://github.com/mrc-ide/mcstate/blob/da9f79e4b5dd421fd2e26b8b3d55c78735a29c27/vignettes/continuous.Rmd#L153
# vignettes/continuous.Rmd
mcmc_sample <- mcstate::pmcmc_sample(pmcmc_tuned_run, 1e3)#1e3

# We can visually check the chains have converged and assess the accuracy of our parameter estimates using a sample from the posterior distribution.
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), bty = "n")
plot(pmcmc_tuned_run$pars[, "beta_0"], type = "l", xlab = "Iteration",
     ylab = "beta_0")
plot(pmcmc_tuned_run$pars[, "beta_1"], type = "l", xlab = "Iteration",
     ylab = "beta_1")
hist(mcmc_sample$pars[, "beta_0"], main = "beta 0", xlab = "beta_0",
     freq = FALSE)
abline(v = pars$beta_0, lty = 2, col = "darkred", lwd = 3)
hist(mcmc_sample$pars[, "beta_1"], main = "beta 1", xlab = "beta_1",
     freq = FALSE)
abline(v = pars$beta_1, lty = 2, col = "darkred", lwd = 3)
legend("topright", legend = "Initial Value", col = "darkred", lty = 2, bty = "n")


par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), bty = "n")
plot(pmcmc_tuned_run$pars[, "wane"], type = "l", xlab = "Iteration",
     ylab = "wane")
plot(pmcmc_tuned_run$pars[, "log_delta"], type = "l", xlab = "Iteration",
     ylab = "log_delta")
hist(mcmc_sample$pars[, "wane"], main = "wane", xlab = "wane",
     freq = FALSE)
abline(v = pars$wane, lty = 2, col = "darkred", lwd = 3)
hist(mcmc_sample$pars[, "log_delta"], main = "log(delta)", xlab = "log_delta",
     freq = FALSE)
abline(v = pars$log_delta, lty = 2, col = "darkred", lwd = 3)
legend("topright", legend = "Initial Value", col = "darkred", lty = 2, bty = "n")

par(mfrow = c(2, 1), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), bty = "n")
plot(pmcmc_tuned_run$pars[, "sigma_2"], type = "l", xlab = "Iteration",
     ylab = "sigma_2")
hist(mcmc_sample$pars[, "sigma_2"], main = "sigma_2", xlab = "sigma_2",
     freq = FALSE)
abline(v = pars$sigma_2, lty = 2, col = "darkred", lwd = 3)
legend("topright", legend = "Initial Value", col = "darkred", lty = 2, bty = "n")
par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), bty = "n")

# We can compare our fitted trajectories to the prevalence data by sampling from the our comparison distribution to obtain an estimate of the number of positive tests under our model.
state <- mcmc_sample$trajectories$state

# Prevalence basically Disease/N
model_prev <-  t(y[4, , -1] / (y[2, , -1] + y[3, , -1] + y[4, , -1] + y[5, , -1]))
modelled_positives <- apply(model_prev, 2, rbinom, n = nrow(incidence),
                            size = incidence$cases) # size should be data$tested (but we don't have tested data (only positive data))

# par(mfrow = c(1, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), bty = "n")
# times <- c(0, incidence$day)
# matplot(times, t(state[2, , ]), type = "l", lty = 1, col = cols["S"],
#         ylim = c(0, max(state)), xlab = "Day", ylab = "Number of individuals")
# matlines(times, t(state[5, , ]), lty = 1, col = cols["R"])
# matlines(times, t(state[4, , ]), lty = 1, col = cols["I"])
# legend("left", lwd = 1, col = cols[-4], legend = names(cols)[-4], bty = "n")

par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), bty = "n")
matplot(incidence$day, modelled_positives, type = "l", lty = 1, col = grey(0.7, 0.2),
        xlab = "Day", ylab = "Number of positive tests")
points(incidence$day, incidence$cases, pch = 20)

# coda::acfplot(mcmc1[, "time_shift"], main = "Autocorrelation time shift")
coda::acfplot(mcmc2[, "beta_0"], main = "Autocorrelation beta0")
coda::acfplot(mcmc2[, "beta_1"], main = "Autocorrelation beta1")
coda::acfplot(mcmc2[, "wane"], main = "Autocorrelation wane")
coda::acfplot(mcmc2[, "log_delta"], main = "Autocorrelation log(delta)")
coda::acfplot(mcmc2[, "sigma_2"], main = "Autocorrelation sigma2")

## 4. Running predictions
