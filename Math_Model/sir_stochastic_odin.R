if (!require(odin, quietly=T)){
  install.packages("odin")
  library(odin)
}

# 28.04.2024
# error occurs when I tried to run the code for both on Win (computer library) & Lin (Mint)
# dde package is required
if (!require(dde, quietly=T)){
  install.packages("dde")
  library(dde)
}

gen_sir <- odin::odin({
  dt <- user()
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
  # 1. PARAMETERS ################################################################
  S_ini <- user()
  I_ini <- user()
  beta <- user()
  sigma <- user()
  
  # 2. INITIAL VALUES ############################################################
  initial(S) <- S_ini
  initial(I) <- I_ini
  initial(R) <- 0
  
  # 3. UPDATES ###################################################################
  N <- S + I + R
  lambda <- beta*I/N
  
  # Individual probabilities of transition
  p_SI <- 1- exp(-lambda * dt) # assume dt equals to 1
  p_IR <- 1- exp(-sigma * dt) # assume dt equals to 1
  
  # Draws for numbers changing between compartments
  n_SI <- rbinom(S, p_SI)
  n_IR <- rbinom(I, p_IR)
  
  # The transitions
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
})

# sir_model$state()
# all_date <- data.frame(date = seq(1:4745)) # is number of day of observations

pars <- list(dt = 1,
  S_ini = 1e5,
  I_ini = 10,
  beta = 0.2,
  sigma = 0.1
  # DOI = 15.75, # 15.75 days (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
)

sir_model <- gen_sir$new(user = pars)

set.seed(0)
timesteps <- seq(0, nrow(all_date), by=1)   # time.
# timesteps <- seq(0, 500, by=1)   # trial other timestep value to compare the result
sir_res <- sir_model$run(timesteps)

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")

matplot(sir_res[, 1], sir_res[, (3:ncol(sir_res))], # not choose the first & second column
        xlab = "Time", ylab = "Number of individuals",
        type = "l", col = cols, lty = 1)
legend("topright", lwd = 1, col = cols, legend = c("S", "I", "R"), bty = "n")

# The sir_res output is a matrix, save this to csv:
write.csv(sir_res, file="Output_sir_result.csv", row.names =T)
