kernmat <- function(x, pars, fi){
  
  ### create a mesh grid of size changes
  Size_c = t(matrix(x, nrow=length(x), ncol=length(x))) #- matrix of sizes at t  ## sizes same per column
  
  ### AGE selectivity eqn A.1.29 [appendix]
  #S = 1 /1 + e(-log(19)(L-B1) / B2)
  ## B1: len @50% female B2: diff btw len 95%-50% length at time t
  selc <- 1/ (1 + exp(-log(19) * (Size_c - pars$selc.50 ) / (pars$selc.100 - pars$selc.50)))
  
  #### Survival part of kernel
  #m=0 #check with mort. =0 
  m= pars$M + fi * selc # natural and fishing
  
  ### convert to probability
  timestep = 1 # 1yr
  pm <- exp(-m * timestep) # prob. survival
  # increases with size
  pm[pm<0] = 0 # make sure no negatives 
  
  ### create a mesh grid of size changes
  Size_r = t(Size_c) #- matrix of sizes at t+1 ## sizes same per row
  
  ### Growth part of kernel 
  nLinfs = 1000 # how many different values of Linf to simulate
  ### create Matrices of simulated sizes
  Linfs <- rnorm(n = nLinfs, mean = pars$Linf, (pars$Linf.cv *pars$Linf)) # vector of random dist. of Linfs- max lengths
  Linf_sims <- matrix(Linfs, nrow = nLinfs, ncol = length(x)) # expand into a matrix so there is a corresponding value of Linf for each possible value in the length vector x
  # ^ same values for each column to run against different size distr.
  Size_chge <- t(matrix(x, ncol = nLinfs, nrow = length(x))) # expand x into a matrix so there is a value for each value of Linfs.mat ## size changes per column and goes up to mesh max size 
  
  ### Create range of growth rates
  VBlengths <-  Linf_sims - (Linf_sims - Size_chge) * exp(-pars$K * timestep) # VonBert -- use those two matrices to get the range of possible growth rates, as a function of X
  
  VBmean <- colMeans(VBlengths) # Take the mean across all of the different trajectories for each value of x
  Leng_mean <- t(matrix(VBmean, nrow = length(x), ncol = length(x))) # expand into a matrix with a corresponding value for each value of Y (the size at time t+1)
  
  VBsd <- VBmean * (pars$Linf.cv)
  Leng_sd <- t(matrix(VBsd, nrow = length(x), ncol = length(x))) 
  
  ### convert to probability 
  pg <- dnorm(Size_r, mean = Leng_mean, sd = Leng_sd) # use dnorm to get the distribution of growth rates 
  pg[pg<0] = 0 # make sure no negatives
  
  ### Run the kernel as prob. death * growth
  K = pm * pg 
  return(K)
}

