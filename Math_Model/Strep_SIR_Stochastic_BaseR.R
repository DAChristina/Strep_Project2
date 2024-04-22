rm(list=ls())

# 1. PARAMETERS & INITIAL VALUES ###############################################
N = 1000
S_int = N - 1
I_int = 1
R_int = 0

R0 = 2 # Assumes R0 equals to 2 (have to check the EpiEstim output for R(t))
DOI = 15.75 # 15.75 days (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
k = 0.1
time_step = 0.1
time_sequence = seq(0, 100, by = time_step)

int_values = c(S_int, I_int, R_int,0,0)
data = array(dim = c(length(time_sequence), 6))
data = data.frame(data)
names(data)= c('time_step', 'S', 'I', 'R', 'inc', 'rec')

data$time_step = time_sequence
data[1,2:6] = int_values

# 3. UPDATE FUNCTIONS ##########################################################
inc <- function(S, I){
  lambda = (I * R0)/(N * DOI)
  return(rbinom(n = 1,prob = lambda * time_step, size = S))   
}

rec <- function(I){
  return(rbinom(n = 1,prob = time_step/DOI, size = I))   
}

model <- function(iterated){
  
  for (i in 2:dim(data)[1]){
    recovery = rec(data$I[(i-1)])
    
    # the new code for incidence ------------------
    if(recovery>0)
      incidence=sum(rnbinom(n=recovery, mu = R0*data$S[i-1]/N, size = k))  #size is the overdispersion parameter k, smaller k -> greater variance
    else
      incidence=0
    incidence = min(incidence,data$S[i-1]) #a check to ensure there are not more incident infections than susceptible individuals
    # ---------------------------------------------
    
    data$S[i] = data$S[(i-1)] - incidence
    
    data$I[i] = data$I[(i-1)] + incidence - recovery
    
    data$R[i] = data$R[(i-1)] + recovery
    
    data$inc[i] = incidence
    data$rec[i] = recovery
  }
  if(iterated) 
    return(data[dim(data)[1],])
  else
    return(data)
}

one_run<-model(iterated = FALSE)
plot(y=one_run$I, x = one_run$time_step, type='l', col="red", ylim=c(0,1000))
lines(y=one_run$S, x = one_run$time_step, col="blue")


iterations <- 100
results = array(dim=c(iterations, 6))
results = data.frame(results)
names(results )= c('Iteration', 'S', 'I', 'R', 'inc', 'rec')

for ( j in 1:iterations){
  results[j,] <- model(iterated = TRUE)
}

hist(results$R)

# plot(x = data$time_step, y = data$inc, type = 'l')
# plot(x = data$time_step, y = data$E, type = 'l')
# plot(x = data$time_step, y = data$I, type = 'l')
