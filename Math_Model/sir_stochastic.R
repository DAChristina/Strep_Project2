

# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
dt <- user()
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
lambda <- (R0*I)/(N*DOI)

# Individual probabilities of transition
p_SI <- 1- exp(-lambda * dt)
p_IR <- 1- exp(-sigma * dt)

# Draws for numbers changing between compartments
# temp_n_SI <- rbinom(S, p_SI) # contrasting to only rbinom(S, p_SI), need correction
n_IR <- rbinom(I, p_IR) # contrasting to only rbinom(I, p_IR), need correction
# size is the overdispersion parameter k, smaller k -> greater variance
# a check to ensure there are not more incident infections than susceptible individuals
# temp_n_SI <- rnbinom(R, R0*S/N)
# n_SI <- if (n_IR > 0) min(temp_n_SI, S)  else 0 # Check if the result < 0, throw 0

n_SI <- if (n_IR > 0) min(rpois(R), S)  else 0


# Error: Expected 2 arguments in rnbinom call, but recieved 3
# Error: Argument to sum must be a symbol or indexed array
# sum(rnbinom(n=recovery, mu = R0*data$S[i-1]/N, size = k))

update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
