# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
freq <- user(1)
dt <- 1/freq
initial(time) <- 0

# 1. PARAMETERS ################################################################
S_ini <- user(0) # required in mcState
A_ini <- user(0) # required in mcState
D_ini <- user(0) # required in mcState
time_shift <- user(0)
beta_0 <- user(0)
beta_1 <- user(0)
vacc <- user(0.9*0.862) # FIXED vaccination coverage * efficacy for infants
log_delta <- user(0) # required in mcState
sigma_1 <- user(1/15.75) # FIXED per-day, carriage duration (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
sigma_2 <- user(1) # FIXED assumed as acute phase
mu_0 <- user(0) # background mortality, assumed as closed system
mu_1 <- user(192/(4064*4745)) # FIXED disease-associated mortality; ratio 192/4064 in 4745 days
# p <- user(0.1) # Reduction in susceptibility post-infection (based on Lilith's)
pi <- user(3.141593) # FIXED

# 2. INITIAL VALUES ############################################################
initial(S) <- S_ini
initial(A) <- A_ini
initial(D) <- D_ini
initial(R) <- 0
initial(n_AD_daily) <- 0
initial(n_AD_cumul) <- 0

# 3. UPDATES ###################################################################
N <- S + A + D + R
beta_temporary <- beta_0*(1+beta_1*sin(2*pi*(time_shift+time)/365))
# Infant vaccination coverage occurs when PCV13 introduced in April 2010 (day 2648 from 01.01.2003)
# https://fingertips.phe.org.uk/search/vaccination#page/4/gid/1/pat/159/par/K02000001/ati/15/are/E92000001/iid/30306/age/30/sex/4/cat/-1/ctp/-1/yrr/1/cid/4/tbm/1/page-options/tre-do-0
beta <- if (time >= 2648) beta_temporary*(1-vacc) else beta_temporary
lambda <- beta*(A+D)/N # infectious state from Asymtomatic & Diseased individuals
delta <- 10^(log_delta)

# Individual probabilities of transition
p_SA <- 1- exp(-lambda * dt)

p_Asym <- 1- exp(-(delta+sigma_1) * dt)
p_AD <- 1- exp(-(delta/(delta+sigma_1) * dt))
# p_AR <- 1- exp(-(sigma_1/(delta+sigma_1) * dt))

p_Dis <- 1- exp(-(sigma_2+mu_0+mu_1) * dt)
p_DR <- 1- exp(-(sigma_2/(sigma_2+mu_0+mu_1)) * dt)
# p_Dd <- 1- exp(-((mu_0+mu_1)/(sigma_2+mu_0+mu_1)) * dt) # disease-associated mortality

# Draws for numbers changing between compartments
n_SA <- rbinom(S, p_SA)
n_Asym <- rbinom(A, p_Asym) # n_Asym <- n_AD + n_AR cause cyclic dependency error
n_AD <- rbinom(n_Asym, p_AD)
n_AR <- n_Asym - n_AD # unless error occurs: Error: 1 particles reported errors. Invalid call to binomial with n = -1, p = 0.0846896, q = 0.91531
n_Dis <- rbinom(D, p_Dis)
n_DR <- rbinom(n_Dis, p_DR)
n_Dd <- n_Dis - n_DR  # unless error occurs: Error: 1 particles reported errors. Invalid call to binomial with n = -1, p = 0.0846896, q = 0.91531

# The transitions
update(time) <- (step + 1) * dt
update(S) <- S - n_SA
update(A) <- A + n_SA - (n_AD + n_AR)
update(D) <- D + n_AD - (n_DR + n_Dd)
update(R) <- R + n_AR + n_DR
# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model
update(n_AD_daily) <- n_AD #if (step %% freq == 0) n_AD else n_AD_daily + n_AD
update(n_AD_cumul) <- n_AD_cumul + n_AD # no interest in asymptomatic cases that've recovered
