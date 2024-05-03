

# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
dt <- user(1) # required in mcState
initial(time) <- 0

# 1. PARAMETERS ################################################################
S_ini <- user(1e5) # required in mcState
I_ini <- user(10) # required in mcState
beta <- user()
sigma <- user()

# 2. INITIAL VALUES ############################################################
initial(S) <- S_ini
initial(I) <- I_ini
initial(R) <- 0
initial(n_SI_daily) <- rbinom(S_ini, p_SI)
initial(n_SI_cumul) <- rbinom(S_ini, p_SI)

# 3. UPDATES ###################################################################
N <- S + I + R
lambda <- beta*I/N

# Individual probabilities of transition
p_SI <- 1- exp(-lambda * dt)
p_IR <- 1- exp(-sigma * dt)

# Draws for numbers changing between compartments
n_SI <- rbinom(S, p_SI)
n_IR <- rbinom(I, p_IR)

# The transitions
update(time) <- (step + 1) * dt
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
# that "little trick" that previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
update(n_SI_daily) <- n_SI
update(n_SI_cumul) <- n_SI + I + n_IR + R
