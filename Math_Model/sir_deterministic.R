if (!require(odin, quietly=T)){
  install.packages("odin")
  library(odin)
}

gen_sir <- odin::odin({
  # dt <- user()
  # initial(time) <- 0
  # update(time) <- (step + 1) * dt
  
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
  p_SI <- 1- exp(-lambda) # assume dt equals to 1
  p_IR <- 1- exp(-sigma) # assume dt equals to 1
  
  # Draws for numbers changing between compartments
  n_SI <- rbinom(S, p_SI)
  n_IR <- rbinom(I, p_IR)
  
  # The transitions
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
})

# sir_model$state()
all_date <- data.frame(date = seq(1:4745)) # is number of day of observations

pars <- list(#dt = 1,
             S_ini = 1e5,
             I_ini = 1,
             beta = 0.2,
             sigma = 0.1
             # DOI = 15.75, # 15.75 days (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
)

sir_model <- gen_sir$new(user = pars)

set.seed(0)
timesteps <- seq(0, nrow(all_date), by=1)   # time.
sir_res <- sir_model$run(0:nrow(all_date))

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")

matplot(x_res[, 1], x_res[, -1], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = cols, lty = 1)
legend("topleft", lwd = 1, col = sir_col, legend = c("S", "I", "R"), bty = "n")