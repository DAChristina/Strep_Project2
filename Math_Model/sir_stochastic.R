# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
dt <- user(1) # required in mcState
initial(time) <- 0

# 1. PARAMETERS ################################################################
S_ini <- user(6e7) # required in mcState
A_ini <- user(100) # required in mcState
D_ini <- user(0) # required in mcState
beta_0 <- user(0.01)
beta_1 <- user(0.9)
vacc <- user(0.9*0.862) # vaccination coverage * efficacy for infants
delta <- user(1/2000) # required in mcState
qu <- user(0.0002)
sigma <- user(1/15.75) # assumed as acute phase, fixed per-day, carriage duration (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
mu_0 <- user(0) # background mortality, assumed as closed system
mu_1 <- user(192/(4064*4745)) # disease-associated mortality; ratio 192/4064 in 4745 days
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
beta_temporary <- beta_0*(1+beta_1*sin(2*pi*(time_shift+time)/365))
# Infant vaccination coverage occurs when PCV13 introduced in April 2010 (day 2648 from 01.01.2003)
# https://fingertips.phe.org.uk/search/vaccination#page/4/gid/1/pat/159/par/K02000001/ati/15/are/E92000001/iid/30306/age/30/sex/4/cat/-1/ctp/-1/yrr/1/cid/4/tbm/1/page-options/tre-do-0
beta <- if (time >= 2648) beta_temporary*(1-vacc) else beta_temporary
lambda <- beta*(A+D)/N # infectious state from Asymtomatic & Diseased individuals

# Individual probabilities of transition
p_SA <- 1- exp(-lambda * dt)
p_Asym <- 1 - exp(-delta * dt)
p_AD <- 1- exp(-(qu * dt)) # cumulative prob of (delta*qu + (delta*(1-qu))) = delta, delta*qu/delta = qu
p_AR <- 1- exp(-((1-qu) * dt)) # cumulative prob of (delta*qu + (delta*(1-qu))) = delta, delta*(1-qu)/delta = (1-qu)
p_DR <- 1- exp(-sigma * dt)
p_Dd <- 1- exp(-(mu_0+mu_1) * dt) # disease-associated mortality

# Draws for numbers changing between compartments
n_SA <- rbinom(S, p_SA)
n_Asym <- rbinom(A, p_Asym) # n_Asym <- n_AD + n_AR
n_AD <- rbinom(n_Asym, p_AD)
n_AR <- rbinom(n_Asym, p_AR)
n_DR <- rbinom(D, p_DR)
n_Dd <- rbinom(D, p_Dd)

# The transitions
update(time) <- (step + 1) * dt
update(S) <- S - n_SA
update(A) <- A + n_SA - (n_AD + n_AR)
update(D) <- D + n_AD - n_DR
update(R) <- R + n_AR + n_DR - n_Dd
# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
update(n_AD_daily) <- n_AD
update(n_AD_cumul) <- n_AD + D + (R - n_AR) # no interest in asymptomatic cases that've recovered
