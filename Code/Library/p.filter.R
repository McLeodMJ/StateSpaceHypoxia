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

p.filter <- function(dat, hypox_a, hypox_b, fi, cv_q, sigma_p, scale, rec1, rec2, rec3, rec4, rec5, rec6, rec7, rec8, rec9, rec10, rec11, rec12, rec13, rec14, rec15, rec16, rec17, rec18, rec19, rec20){

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
    Avg.likel <- rep(NA, time) #average of the ftmp
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
  det <- inv_logit(hypox * hypox_a + hypox_b)
  
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
    pois_mle <- data.frame(lambda_vals = Nf[,t,q] * dx * WClen) %>%
      rowwise() %>%
      mutate(log_likelihood = pois.likel(y = Nact[,t], lambda = lambda_vals)) %>%
      ungroup() # avoids issues with dpois() over vectorized data
    
    likel[,q] <- pois_mle$log_likelihood
  } # end of Q loop
    likel[which(!is.finite(likel))] <- log(1e-320) # smallest non-infinite value
    
      
      ## factoring in detection parameter to likelihood data
       if(sum(Nact[,t]) == 0){ 
         # absent
        ftmp[t,] <- colSums(likel,na.rm=TRUE) * c(1 - det[t]) / c(sum(det) * time)
        }else{
          #present
         ftmp[t,] <- colSums(likel,na.rm=TRUE) * det[t] /c(sum(det) * time)
        }
      # shoud detection be going through a likelihood function of its own or just multiplying by the poislike is fine??
      Avg.likel[t] = mean(ftmp[t,])
      Avg.likel[which(!is.finite(Avg.likel))] <- log(1e-320) # smallest non-infinite value
      
   ## Step 4: Resample    
    Wgt = cumsum(ftmp[t,]/ sum(ftmp[t,]))
    Rnd = runif(Q)
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
# dat = Pop.eq$Pop.matrix
# Nact = Pop.sim
# fish = "Dsole"
# hypox = runif(20)*3
# hypox_p = 4  
# fi = Dsole_f
# mesh = 200
# Q = 100
# time = time
# cv = 0.01
# rec_var = rep(0.1, time)
# cv_q = 0.1
# sigma_p = 0.1
# scale = 0.00001
# 