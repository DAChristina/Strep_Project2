# A nice Intro & some examples:
# https://github.com/mrc-ide/odin-dust-tutorial/
# https://mrc-ide.github.io/odin.dust/articles/sir_models.html
# https://mrc-ide.github.io/sircovid/

if (!require(odin.dust, quietly=T)){
  install.packages("odin.dust",
                   repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
  
  library(odin.dust)
}

# This will be saved as 'sir_stochastic.R' and will be run by using library(odin.dust)
# Set wd to the saved file or sir_stochastic.R:
rm(list=ls())

wd = "C:/Users/dac23/Documents/Downloads" # DIDE
wd = "C:/Users/dac23/Downloads" # library computers
wd = "/home/ron/Downloads" # personal OSs
setwd(wd)

gen_sir <- odin.dust::odin_dust("sir_stochastic.R")

# Running the SIR model with dust
pars <- list(dt = 1,
             S_ini = 6e7, # England's pop size is roughly 67,000,000
             A_ini = 100,
             D_ini = 0,
             beta_0 = 0.25, # based on mcmc run, prior = (1/15.75)*1.2,
             beta_1 = 0.8,
             delta = (1/60), # incubation period, trial random number
             qu = 0.23, # When clinical symptoms occur, (1-qu) = no clinical symptoms detected
             sigma = (1/15.75) # carriage duration of diseased = 15.75 days (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
) # Serotype 1 is categorised to have the lowest carriage duration

sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = 1L,
                         n_threads = 4L,
                         seed = 1L)

sir_model$state()

# update_state is required "every single time" to run & produce matrix output (don't know why)
sir_model$update_state(pars = pars,
                       time = 0) # make sure time is 0

# all_date <- incidence$day
# all_date <- data.frame(col = integer(4745))
n_times <- 4745 # 4745 or similar to the number of date range (of the provided data), or try 500 for trial
n_particles <- 10
x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

for (t in seq_len(n_times)) {
  x[ , , t] <- sir_model$run(t)
}
time <- x[1, 1, ] # because in the position of [1, 1, ] is time
x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
library(tidyverse)
glimpse(x)

par(mar = c(5.1, 5.1, 0.5, 0.5), mgp = c(3.5, 1, 0), las = 1)
cols <- c(S = "#8c8cd9", A = "darkred", D = "#cc0099", R = "#999966", n_AD_daily = "orange", n_AD_cumul = "green")
# matplot(time, t(x[1, , ]), type = "l",
#         xlab = "Time", ylab = "Number of individuals",
#         col = cols[["S"]], lty = 1, ylim = range(x))
matplot(time, t(x[2, , ]), type = "l",
                xlab = "Time", ylab = "Number of individuals",
                col = cols[["A"]], lty = 1)#, ylim = max(x[2,,]))
# matlines(time, t(x[2, , ]), col = cols[["A"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["D"]], lty = 1)
# matlines(time, t(x[4, , ]), col = cols[["R"]], lty = 1)
# matlines(time, t(x[5, , ]), col = cols[["n_AD_daily"]], lty = 1)
# matlines(time, t(x[6, , ]), col = cols[["n_AD_cumul"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
max(x[5,,]) # Check max n_AD_daily
max(x[3,,]) # Check max D
pars$beta_0/(pars$delta) + (pars$qu*pars$delta)/(pars$delta*pars$sigma) # print R0

write.csv(x, file="Output_sir_result.csv", row.names =T)
