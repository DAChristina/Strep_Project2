

# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
dt <- user(1) # required in mcState
initial(time) <- 0

# 1. PARAMETERS ################################################################
S_ini <- user(1e5) # required in mcState
A_ini <- user(0) # required in mcState
I_ini <- user(10) # required in mcState
beta <- user()
delta <- user(0.2) # required in mcState
sigma <- user(1/15.75) # fixed per-day, carriage duration (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
pi <- user(3.141593)

# 2. INITIAL VALUES ############################################################
initial(S) <- S_ini
initial(A) <- A_ini
initial(I) <- I_ini
initial(R) <- 0
initial(n_AI_daily) <- rbinom(A_ini, p_AI)
initial(n_AI_cumul) <- rbinom(A_ini, p_AI)

# 3. UPDATES ###################################################################
N <- S + A + I + R
lambda <- beta*I/N

# Individual probabilities of transition
p_SA <- 1- exp(-lambda * dt)
p_AI <- 1- exp(-delta * dt)
p_IR <- 1- exp(-sigma * dt)

# Draws for numbers changing between compartments
n_SA <- rbinom(S, p_SA)
n_AI <- rbinom(A, p_AI)
n_IR <- rbinom(I, p_IR)

# The transitions
update(time) <- (step + 1) * dt
update(S) <- S - n_SA
update(A) <- A + n_SA - n_AI
update(I) <- I + n_AI - n_IR
update(R) <- R + n_IR
# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
update(n_AI_daily) <- n_AI
update(n_AI_cumul) <- n_AI + I + n_IR + R
