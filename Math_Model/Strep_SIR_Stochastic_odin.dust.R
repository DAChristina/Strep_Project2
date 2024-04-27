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

library(odin.dust)
gen_sir <- odin.dust::odin_dust("sir_stochastic.R")

# Running the SIR model with dust
sir_model <- gen_sir$new(pars = list(dt = 1,
                                     S_ini = 1e5,
                                     I_ini = 1,
                                     # DOI = 15.75, # 15.75 days (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
                                     beta = 0.2,
                                     sigma = 0.1
                                     ),
                         time = nrow(all_date),
                         n_particles = 1L,
                         n_threads = 4L,
                         seed = 1L)

# sir_model$state()

n_times <- nrow(all_date) # day
x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

for (t in seq_len(n_times)) {
  x[ , , t] <- sir_model$run(t)
}
time <- x[1, 1, ]
x <- x[-1, , ]

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = range(x))
matlines(time, t(x[2, , ]), col = cols[["I"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["R"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
