####################################
# particle filter function
#' dat: simulated IPM data --> equilibrium
#' Nact: actual data from the WCGBTS [if checking work - then this is 2nd run of IPM]
#' fish: "Dsole", "Lingcod", "Yrock", "Grock" to call param for species of interest
#' fi: fishing rate - pars$f generally fits to call from param() fxn
#' hypox: hypoxia parameter - need to adjust for MCMC
#' selc: includes WCGBTS selectivity. if do not want to include selc = NA
#' mesh: mesh size - defines nrows in matrix
#' Q: no. of particles 
#' time: running iterations to this time
#' cv: variation for recruitment variation
#' cv_q: random variation introduced for each particle 
#' sigma_p: variation for prcoess error

#####################################
### PraCTICE IF WORKS

#load("../Data/Hist_spp.Rdata")
#load("./Results/Expected_Fishing_rate.Rdata")
#source("./Library/params.R")
#source("./Library/kernmat.R")
#source("./Library/fecmat.R")

#######################################

p.filter <- function(dat, Nact, fish, hypox, fi, selc, mesh, Q, time, cv, cv_q, sigma_p){
  
  # Particle filter: (following Knape & deValpine 2012)
  N <- matrix(NA, nrow = nrow(dat), ncol = time)
  
  # Generate Q particles (independent simulations of N):
  Nf = array(rep(NA, nrow(N) * time * Q), c(nrow(N), time, Q)) #pop.dist. x time x particles
  Nf[,1,] <- dat[,time] #saves last column of equil. data
  
  # Simulate random  variation for first particle
  Nf[,1,] <- Nf[,1,] + rnorm(nrow(N), 0, cv_q * dat[,time])
  Nf[Nf<0] = 0 # check for all values greater than zero 
  
  # Step 1: Initialize the model 
    # resample for accuracy
    ftmp <- matrix(NA, nrow= time, ncol = Q)
    ftmp[1,] <- rep(1, Q)
    Wgt = cumsum(ftmp[1,]/ sum(ftmp[1,]))  
    Rnd = runif(Q)
    Wgt = matrix(Wgt, nrow=Q, ncol=Q) # same across row?
    Rnd =  t(matrix(Rnd, nrow=Q, ncol=Q)) # same across column?
    Pass <- matrix(Rnd < Wgt, nrow=Q, ncol= Q)
    Ind = Q - colSums(Pass) + 1
    
  Nf[,1,] <- Nf[ ,1, Ind] # replace with resampled values -- based on Ind variable that pulls the best representable samples 
  
  Ind2 = sample(1:Q, 1) # index so must be less than 100
  N[,1] = Nf[,1,Ind2] # pick one randomly to be *the* distribution to carry forward to the next step
  
  
  ## Step 2: run model through time T 
    # Create the kernel:
    pars <- params(fish)
    
    meshmax = pars$Linf * 2
    x <- seq(from = 1, to = meshmax, length.out = mesh)
    dx = diff(x)[1]
    
  ## size distribution of recruits
  Rvec <- dnorm(x, pars$Rlen, pars$Rlen.sd)  
  
  K <- kernmat(x, pars, Dsole_f,  1)
  Fe <- fecmat(x, pars)
  
  for (t in 2:time){
    
    # Advance the model, one particle at a time
    for(q in 1:Q){
      Nrand = rnorm(nrow(N), 0 , sigma_p) #processs error
      
      # integration
      Nf[,t,q] <- K %*% Nf[,t-1,q]  * dx  # midpoint rule integration
      E <- sum(Fe * Nf[,t-1,q]) * dx 
      
      # Beverton-holt density-depedence
      Recruits <- as.vector((4 * pars$steep * exp(pars$R0) * E) / ((pars$S0 * (1 - pars$steep)) + (E * (5 * pars$steep - 1))))
      # Add variation in recruitment 
      RR <- exp(rnorm(1, mean = log(Recruits) - ((cv * log(Recruits))^2)/2, sd= cv * log(Recruits) ))# change cv to 0 for no variation
      
      Nf[,t, q] = Nf[,t, q] + (Rvec * RR) + Nrand  # move the model forward for each particle 
    } # end of Q loop
    Nf[Nf < 0] = 0
    
    # Need to look into WCGBTS Selectivity
    if (length(selc) == 1){
       likel <- dpois(matrix(rep(Nact[,t], Q), ncol=Q), Nf[,t,] * dx, log= T) 
      likel[which(!is.finite(likel))] <- 0
      ftmp[t,] = colSums(likel) #log-likelihood vector
    } else {
      #OKlen <- selc * dat
      likel <- dpois(matrix(rep(Nact[,t], Q), ncol=Q), Nf[,t,] * dx * rep(WClen[,Q]), log= T) 
      likel[which(!is.finite(likel))] <- 0
      ftmp[t,] = colSums(likel) #log-likelihood vector
      #ftmp(t,:) = sum(log(max(realmin,poisspdf(repmat(Nact(:,t),[1,1,Q]),(NT(t-length(Tpre))*dy*Nf(:,t,:)).*repmat(OKlen(:),[1,1,Q])))));
      # oklen = filter out sizes wouldnt see in the data (btw 0-1: with full selection at 1)
    }
    
    # create space for false positives
    enc.2 <- rbinom(time, 1, hypox/10) # hypox rate 
    ftmp.temp <- ftmp #replicate 
    for(i in 1:time){
      if(enc.2[i] == 0){
        ftmp.temp[i,] <- rep(-999, Q) # replace zero observ. moments with -999
      }}
    
    
    # if averages over trawl (lognormal) otherwise poisson dpois()
    # repat mat in matrix for replicating actual dat for Q particles
    #' oklen Give selectivty - selectvity of wcgbts 
    #' mean is intergration of NF to get counts b/c working with densities prior - midpt rule 
    
    # Resample    
    Wgt = cumsum(ftmp[t,]/ sum(ftmp[t,]))
    Rnd = runif(Q)
    Wgt = matrix(Wgt, nrow=Q, ncol=Q) # same across row?
    Rnd =  t(matrix(Rnd, nrow=Q, ncol=Q)) # same across column?
    Pass = matrix(Rnd < Wgt, nrow=Q, ncol= Q)
    Ind = Q - colSums(Pass) + 1
    
    
    Nf[,t,] <- Nf[ ,t, Ind] # replace with resampled values
    Ind2 = sample(1:Q, 1)
    N[,t] = Nf[,t,Ind2] # pick one randomly to be *the* distribution to carry forward to the next step
    } # end of t loop
  
  
  
  Pfilter_data <- list (likelihood = ftmp.temp,
                         Pop = N)
  
  return(Pfilter_data)
}
