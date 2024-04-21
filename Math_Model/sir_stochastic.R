
# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

# 1. PARAMETERS ################################################################
N <- S + I + R
I0 <- user()
beta <- user()
sigma <- user()

# 2. INITIAL VALUES ############################################################
initial(S) <- N - I0
initial(I) <- I0
initial(R) <- 0


# 3. UPDATES ###################################################################
p_SI <- 1 - exp(-beta * I / N)
p_IR <- 1 - exp(-sigma)

n_SI <- rnbinom(S, p_SI)
n_IR <- rnbinom(I, p_IR)

update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR