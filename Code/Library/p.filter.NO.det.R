####################################
# particle filter function
#' dat: simulated IPM data --> equilibrium
#' Nact: actual data from the WCGBTS [if checking if p.filter works - then this is 2nd run of IPM]
#' selc: F means we are running simulations and do NOT need the WCGBTS selc- T means we need WCGBTS selc. 
#' fish: "Dsole", "Lingcod", "Yrock", "Grock" to call param for species of interest
#' fi: fishing rate - take fishing rate found from SPR Analysis OR pars$f generally fits to call from param() fxn
#' hypox_a: dep. hypoxia parameter - based on the det. prob. that gives at least 4 zeros when run through rbinom() hypox_a=3
#' hypox_b: intercept hypoxia parameter - if DO = 0 & prob is very low [0.001] logit(0.001/.999) --> hypox_b = -7
#' mesh: mesh size - defines nrows in matrix
#' Q: no. of particles 
#' time: running iterations to this time
#' rec_var: variation for recruitment variation
#' cv_q: random variation introduced for each particle 
#' sigma_p: variation for process error
#' scale: used to reduce pop sizes down to more realistic values per trawl - best scale = 0.00001

#######################################

p.filter.NO.det <- function(dat, fi, cv_q, sigma_p, scale, rec1, rec2, rec3, rec4, rec5, rec6, rec7, rec8, rec9, rec10, rec11, rec12, rec13, rec14, rec15, rec16, rec17, rec18, rec19, rec20){
  
  # return data from list
  Nact = dat$N.act
  Nint = dat$Nint
  hypox = dat$DO 
  fish = dat$fish
  mesh = dat$Mesh
  Q = dat$Q
  time = dat$time
  rec_var = c(rec1, rec2, rec3, rec4, rec5, rec6, rec7, rec8, rec9, rec10, rec11, rec12, rec13, rec14, rec15, rec16, rec17, rec18, rec19, rec20)
  
  
  # format data to be scaled down to more realistic sizes
  Nint <- Nint * scale
  Nact <- round(Nact * scale)
  
  # Particle filter: (following Knape & deValpine 2012)
  N <- matrix(NA, nrow = mesh, ncol = time)
  
  # Generate Q particles (independent simulations of N):
  Nf <- array(rep(NA, mesh * time * Q), c(mesh, time, Q)) #pop.dist. x time x particles
  Nf[,1,] <- Nint #inputs the last column of pop. equil. data
  
  # Simulate random variation for first particle
  part_sims <- rnorm(length(Nf[,1,]), 0, cv_q)
  dim(part_sims) <- dim(Nf[,1,]) 
  Nf[,1,] <- Nf[,1,] + part_sims
  Nf[Nf<0] = 0 # check for all values greater than zero 
  
  
  # Step 1: Initialize the model 
  # resample for accuracy
  Avg.likel <- rep(NA, time-1) #average of the ftmp
  likel <- matrix(NA, nrow=mesh, ncol= Q)
  ftmp <- matrix(NA, nrow= time, ncol = Q)
  ftmp[1,] <- rep(1, Q)
  Wgt = cumsum(ftmp[1,]/ sum(ftmp[1,]))  
  Rnd = runif(Q)
  Wgt = matrix(Wgt, nrow=Q, ncol=Q) # same across row
  Rnd =  t(matrix(Rnd, nrow=Q, ncol=Q)) # same across column
  Pass <- matrix(Rnd < Wgt, nrow=Q, ncol= Q) # logical T or F for if the random value is less than the weighted value
  Ind = Q - colSums(Pass) + 1 #pulls the best representable samples 
  
  Nf[,1,] <- Nf[ ,1, Ind] # replace with resampled values
  
  Ind2 = sample(1:Q, 1) # index so must be less than 100
  N[,1] = Nf[,1,Ind2] # pick one randomly to be *the* distribution to carry forward to the next step
  
  ## Step 2: create the operating model
  
  # calls parameters from SS
  pars <- params(fish)
  
  # integration
  meshmax = pars$Linf * 2
  x <- seq(from = 0, to = meshmax, length.out = mesh)
  dx = diff(x)[1]
  
  ## size distribution of recruits
  Rvec <- dnorm(x, pars$Rlen, pars$Rlen.sd)  
  
  # Create the kernel
  K <- kernmat(x, pars, exp(fi)) # very small...
  Fe <- fecmat(x, pars)
  
  # process error
  Nrand <- rnorm(length(Nf[,1,]), 0, exp(sigma_p)) 
  dim(Nrand) <- dim(Nf[,1,])
  #Nrand = matrix(rnorm(mesh * Q, 0 , sigma_p), mesh, Q) #processs error
  
  # Detection parameter dependent on hypoxia parameter and DO data
  #det <- inv_logit(hypox * hypox_a + hypox_b)
  
  # WCGBTS Selectivity
  WClen <- 1 # does NOT add trawl selc. 
  
  
  ## Step 3: run model through time T 
  Recruits <- NULL
  for (t in 2:time){
    
    # Advance the model, one particle at a time
    for(q in 1:Q){
      # integration
      Nf[,t,q] <- K %*% Nf[,t-1,q]  * dx  # midpoint rule integration
      E <- sum(Fe * Nf[,t-1, q] * dx )
      
      # Beverton-holt density-depedence
      Recruits[t] <- as.vector((4 * pars$steep * exp(pars$R0) * E)/ ((pars$S0 * (1 - pars$steep)) + (E * (5 * pars$steep - 1))) ) 
      # THIS gives a reasonable value but does it make sense to apply this way??? - as.vector((4 * pars$steep * exp(pars$R0*scale)* E)/ ((pars$S0*scale * (1 - pars$steep)) + (E * (5 * pars$steep - 1))) ) 
      Nf[,t, q] = Nf[,t, q] +  exp(rec_var[t] + log(Rvec * Recruits[t])) + Nrand[ ,q] # add var. in recruitment
      Nf[Nf < 0] = 0
      
      
      ### site by year analysis 
      
      ############################### LIKELIHOOD ############################
   # pois_mle <- data.frame(lambda_vals = Nf[,t,q] * dx * WClen) %>%
   #   rowwise() %>%
   #   mutate(log_likelihood = pois.likel(y = Nact[,t], lambda = lambda_vals)) %>%
   #   ungroup() # avoids issues with dpois() over vectorized data
   # likel[,q] <- pois_mle$log_likelihood
     
      likel[,q] <- dpois(x=Nact[,t], lambda = Nf[,t,q] * dx * WClen, log=T)
      
    } # end of Q loop
    
    like.min <- min(likel[is.finite(likel)], na.rm = T) 
    #add ifelse statement for when bounds are reached and entire likel matrix returns NaNs, like.min will return Inf
    if(like.min == Inf){
      like.min = -1e50
    }
    likel[which(!is.finite(likel))] <-  like.min  
    
  # taking the sum of the likel over time but NOT looking at det
      ftmp[t,] <- colSums(likel,na.rm=TRUE)
    
    # shoud detection be going through a likelihood function of its own or just multiplying by the poislike is fine??
    Avg.likel[t-1] = mean(ftmp[t,])
   
    ## Step 4: Resample    
    Wgt = cumsum(na.omit(ftmp[t,])/ sum(ftmp[t,], na.rm =T))
    Rnd = na.omit(runif(Q))
    Wgt = matrix(Wgt, nrow=Q, ncol=Q) # same across row
    Rnd =  t(matrix(Rnd, nrow=Q, ncol=Q)) # same across column
    Pass = matrix(Rnd < Wgt, nrow=Q, ncol= Q)
    Ind = Q - colSums(Pass) + 1
    
    Nf[,t,] <- Nf[ ,t, Ind] # replace with resampled values
    Ind2 = sample(1:Q, 1)
    N[,t] = Nf[,t,Ind2] # pick one randomly to be *the* distribution to carry forward to the next step
  } # end of t loop
  
  Pfilter_data <- list(likelihood = -1 *sum(Avg.likel, na.rm = T),
                       Fit = N,
                       True=Nact,
                       recruits = Recruits)
  
  return(Pfilter_data$likelihood)
}
