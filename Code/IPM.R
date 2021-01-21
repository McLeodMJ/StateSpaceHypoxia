# Greenstriped rockfish 
# Fishing rate is referred to as exploitation rate

source("params.R")

# Here's that kernel function
kernmat <- function(x, pars, timestep){
  
  #' create a mesh grid of size changes
   X = t(matrix(x, nrow=length(x), ncol=length(x))) #- matrix of sizes at t
   Y = t(X) #- matrix of sizes at t+1
  
  # Survival part of kernel
  # m= pars$m
  m=0  #check with mort. =0 - # should be dist. around Linf
  
  
  ### decreases with size class

  # convert to probability
  pm <- exp(-m * timestep) # prob. survival
  ### increases with size
  
  # Growth part of kernel. Do it this way so that we simulate many different growth trajectories in the population, averaged together.
  nLinfs = 1000 # how many different values of Linf to simulate
  Len_sims <- array(rep(NA, nLinfs * length(x) * length(x)), c(nLinfs, length(x), length(x)))
  
  Linfs <- rnorm(n = nLinfs, mean = pars$Linf, sd = pars$Linf.sd) # vector of distribution of Linfs
  Len_sims[,,1] <- matrix(Linfs, nrow = nLinfs, ncol = length(x)) # expand into a matrix so there is a corresponding value of Linf for each possible value in the length vector x
  Len_sims[,1,] <- t(matrix(x, ncol = nLinfs, nrow = length(x))) # expand x into a matrix so there is a value for each value of Linfs.mat

  g.tmp <-  Len_sims[,,1] - (Len_sims[,,1] - Len_sims[,1,]) * exp(-pars$k * timestep) # use those two matrices to get the range of possible growth rates, as a function of X
  g.mean <- colMeans(g.tmp) # Take the mean across all of the different trajectories for each value of x
  Len_sims[1,,] <- t(matrix(g.mean, nrow = length(x), ncol = length(x))) # expand into a matrix with a corresponding value for each value of Y (the size at time t+1)
  pg <- dnorm(Y, mean = Len_sims[1,,], sd = 2 * dx) # use dnorm to get the distribution of growth rates (using an arbitrarily small sd)
  
  
  # make sure no negatives
  pg[pg<0] = 0
  pm[pm<0] = 0 
  
  # you would add fecundity here if it were a closed population...
  k = pm*pg
  # get the diagonal growth pattern
  
  # If adding fecundity too:
  #Rvec = dnorm(Y, pars$rec, pars$rec.sd) # size distribution of recruits
  #Fec = pars$fec.const * X ^pars$fec.exp # you can use values of 1 and 3 for these... 
  #Q = Fec * Rvec # combine the two
  #k = k + Q # add in fecundity to the growth/survival part
  
  
  return(k)
}

#------------------------------------------------------------------
# Build IPM...

# Create parameters:
pars <- params("Grock")
  #list(k = 0.105, Linf = 33.67, Lvar = 0.14/33.67, Rmu = 10, # Lvar is CV of Linf (SD/mean)
             # Rmean = 9.65, Rsd = 0.15, mort= 0.0682) #, Fec_constant = 1, Fec_exponent = 3) # mean & sd of size of new recruits

# rates per year 
# linf cm
# nat. mort - informative prior with life history  (max observed age - average water temp --> life history theory )
# size of a new recruits - back estimate  - size at age 0 after von bertalnffy 
# fishing rate 
# variabilty in recruitment 
# can we estimate the fishing rate getting info from edge while acccounting for hypoxia ?


# IPM integration parameters:
meshsize = 100
meshmin = 0
meshmax = pars$Linf * 2

x <- seq(from = meshmin, to = meshmax, length.out = meshsize) # different size classes
dx = diff(x)[1] #width of the 'rectangle' to do the midpoint rule *cough cough* left rule

# Recruit size distribution
Rvec <- dnorm(x, mean = pars$rec, pars$rec.sd)
# this should integrate to unity (1) - does it?

# Create the kernel
K <- kernmat(x, pars, timestep = 1/36) # in 1/36 of a year time step
plot(K, type='l')

# Initialize the model:
time = 1000
N = matrix(0, nrow = meshsize, ncol = time)
N[,1] <- Rvec # * pars$Rmu # initialize with one pulse of new recruits

# Run the model
for (t in 2:time){
  N[,t] <- K %*% N[,t-1] * dx  +  Rvec # * pars$Rmu # midpoint rule integration
}

# total population size:
sum(N[ ,time] * dx)


plot(x, N[ ,time], type = 'l')
plot(colSums(N) * dx)

