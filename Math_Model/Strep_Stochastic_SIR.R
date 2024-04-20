# A nice Intro & some examples:
# https://github.com/mrc-ide/odin-dust-tutorial/
# https://mrc-ide.github.io/odin.dust/articles/sir_models.html

if (!require(odin.dust, quietly=T)){
  install.packages("odin.dust",
                   repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
  
  library(odin.dust)
}

transition <- odin.dust::odin_dust({
# 1. PARAMETERS ################################################################
  N <- user()
  I0 <- user()
  beta <- user()
  sigma <- user()
  
# 2. INITIAL VALUES ############################################################
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  
  n_SI <- rnbinom(S, 1 - exp(-beta * I / N))
  n_IR <- rnbinom(I, 1 - exp(-sigma))

# 3. UPDATES ###################################################################
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  
})

pars <- list(N = 10e6,
             I0 = 1,
             beta = 1,
             sigma = 1)

mod <- transition$new(user = pars) # changing the cycle by user loops instead of define the cycle_width one-by-one
timesteps <- seq(0, 2000, by=1)   # time.
y <- mod$run(timesteps)
