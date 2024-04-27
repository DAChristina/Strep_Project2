

# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
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
p_SI <- 1- exp(-lambda * dt)
p_IR <- 1- exp(-sigma * dt)

# Draws for numbers changing between compartments
n_SI <- rbinom(S, p_SI)
n_IR <- rbinom(I, p_IR)

# The transitions
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
