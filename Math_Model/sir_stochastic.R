# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
dt <- user(1) # required in mcState
initial(time) <- 0

# 1. PARAMETERS ################################################################
S_ini <- user(6e7) # required in mcState
A_ini <- user(0) # required in mcState
D_ini <- user(10) # required in mcState
beta_0 <- user(0.5)
beta_1 <- user(0.01)
delta <- user(0.2) # required in mcState
qu <- user(1/15.75) # fixed per-day, carriage duration (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
sigma <- user((1/15.75)*3) # assumed as acute phase
pi <- user(3.141593)
time_shift <- user(70)

# 2. INITIAL VALUES ############################################################
initial(S) <- S_ini
initial(A) <- A_ini
initial(D) <- D_ini
initial(R) <- 0
initial(n_AD_daily) <- rbinom(A_ini, p_AD)
initial(n_AD_cumul) <- rbinom(A_ini, p_AD)

# 3. UPDATES ###################################################################
N <- S + A + D + R
beta <- beta_0*(1+beta_1*sin(2*pi*(time_shift+time)/365))
lambda <- beta*D/N

# Individual probabilities of transition
p_SA <- 1- exp(-lambda * dt)
p_Asym <- 1 - exp(-delta * dt)
p_AD <- 1- exp(-((delta*qu/delta) * dt)) # cumulative prob of (delta*qu + (delta*(1-qu))) = delta
p_AR <- 1- exp(-((delta*(1-qu)/delta) * dt)) # cumulative prob of (delta*qu + (delta*(1-qu))) = delta
p_DR <- 1- exp(-sigma * dt)

# Draws for numbers changing between compartments
n_SA <- rbinom(S, p_SA)
n_Asym <- rbinom(A, p_Asym) # n_Asym <- n_AD + n_AR
n_AD <- rbinom(n_Asym, p_AD)
n_AR <- rbinom(n_Asym, p_AR)
n_DR <- rbinom(D, p_DR)

# The transitions
update(time) <- (step + 1) * dt
update(S) <- S - n_SA
update(A) <- A + n_SA - (n_AD + n_AR)
update(D) <- D + n_AD - n_DR
update(R) <- R + n_AR + n_DR
# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
update(n_AD_daily) <- n_AD
update(n_AD_cumul) <- n_AD + D + (R - n_AR) # no interest in asymptomatic cases that've recovered
