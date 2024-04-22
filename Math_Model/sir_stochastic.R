
# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

# 1. PARAMETERS ################################################################
S_ini <- user()
I_ini <- user()
R0 <- user()
DOI <- user()
sigma <- user()
# kappa <- user() # overdispersion parameter, with smaller k -> greater variance

# 2. INITIAL VALUES ############################################################
initial(S) <- S_ini
initial(I) <- I_ini
initial(R) <- 0

# 3. UPDATES ###################################################################
N <- S + I + R
beta <- (R0*I)/(N*DOI)

# Individual probabilities of transition
p_SI <- 1- exp(-beta)
p_IR <- 1- exp(-sigma)

# Draws for numbers changing between compartments
temp_n_SI <- rnbinom(S, p_SI)
temp_n_IR <- rnbinom(I, p_IR)
# Error: Expected 2 arguments in rnbinom call, but recieved 3
# temp_n_SI <- rnbinom(n=n_IR, mu = mu, size = kappa) # (line 30)
# size is the overdispersion parameter k, smaller k -> greater variance
# a check to ensure there are not more incident infections than susceptible individuals

n_SI <- if (temp_n_SI > 0) temp_n_SI else 0 # Check if the result < 0, throw 0
n_IR <- if (temp_n_IR > 0) temp_n_IR else 0 # Check if the result < 0, throw 0

update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
