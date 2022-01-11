####################################
# particle filter function
#' dat: simulated IPM data --> equilibrium
#' Nact: actual data from the WCGBTS [if checking if p.filter works - then this is 2nd run of IPM]
#' selc: F means we are running simulations and do NOT need the WCGBTS selc- T means we need WCGBTS selc. 
#' fish: "Dsole", "Lingcod", "Yrock", "Grock" to call param for species of interest
#' fi: fishing rate - take fishing rate found from SPR Analysis OR pars$f generally fits to call from param() fxn
#' hypox: hypoxia parameter - need to adjust for MCMC
#' mesh: mesh size - defines nrows in matrix
#' Q: no. of particles 
#' time: running iterations to this time
#' rec_var: variation for recruitment variation
#' cv_q: random variation introduced for each particle 
#' sigma_p: variation for process error
#' scale: used to reduce pop sizes down to more realistic values per trawl - best scale = 0.00001

#######################################

p.filter.WCGBTS <- function(dat, nact, hypox, hypox_p, fi, cv_q, sigma_p, scale, rec1, rec2, rec3, rec4, rec5){
  # return data from list
  Nint = dat$Nint #inputs the last column of steady state pop.
  fish = dat$fish
  mesh = dat$Mesh
  Q = dat$Q
  rec_var <- list("2010" = rec1, 
                     "2011" = rec2,
                     "2012" = rec3, 
                     "2013" = rec4, 
                     "2014" = rec5)

  
  # format data to be scaled down to more realistic sizes
  Nint <- Nint * scale
  #Nact <- round(Nact * scale)
  
  ### Step 1:  site by year analysis
  Avg.likel = NULL
  # empty list for each yrs simulated pop using p.filter
  Pop_list <- vector(mode = "list", length = length(nact))
  names(Pop_list) <- names(nact)
  
  ## additional for loop to run over each year & trawl sites as columns
  for(y in length(nact)){
  y=1
    hyp <- hypox[[y]]$mean_hyp
    # Detection parameter dependent on hypoxia parameter and DO data
    det <- inv_logit(hyp * hypox_p)
    
    #population by year
    Nact <- nact[[y]]
    Nact <- round(Nact * scale) # if we want to scale it down 
    # do we want to return each year or take the avg of the (-sum LLs?)
    
    # time must be formatted by matrices of each year (ncol - trawl sites)
    time = ncol(Nact)
  
  # Step 2: Particle filter: (following Knape & deValpine 2012)
  N <- matrix(NA, nrow = mesh, ncol = time)
  
  # Generate Q particles (independent simulations of N):
  Nf <- array(rep(NA, mesh * time * Q), c(mesh, time, Q)) #pop.dist. x time x particles
  Nf[,1,] <- Nint #inputs the last column of pop. equil. data
  
  # Simulate random variation for first particle
  part_sims <- rnorm(length(Nf[,1,]), 0, cv_q)
  dim(part_sims) <- dim(Nf[,1,]) 
  Nf[,1,] <- Nf[,1,] + part_sims
  Nf[Nf<0] = 0 # check for all values greater than zero 
  
  
  # Step 3: Initialize the model 
  # resample for accuracy
  avg.likel <- rep(NA, time) #average of the ftmp
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
  
  
  ## Step 4: create the operating model
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
  
  # WCGBTS Selectivity
   WClen <- pars$selc 
  
   
  ## Step 5: run model through time T 
  Recruits <- NULL
 for (t in 2:time){
    
    # step 5: Advance the model, one particle at a time
    for(q in 1:Q){
      # integration
      Nf[,t,q] <- K %*% Nf[,t-1,q]  * dx  # midpoint rule integration
      E <- sum(Fe * Nf[,t-1, q] * dx )
      
      # Beverton-holt density-depedence
      Recruits[t] <- as.vector((4 * pars$steep * exp(pars$R0) * E)/ ((pars$S0 * (1 - pars$steep)) + (E * (5 * pars$steep - 1))) ) 
      # THIS gives a reasonable value but does it make sense to apply this way??? - as.vector((4 * pars$steep * exp(pars$R0*scale)* E)/ ((pars$S0*scale * (1 - pars$steep)) + (E * (5 * pars$steep - 1))) ) 
      Nf[,t, q] = Nf[,t, q] +  exp(rec_var[[y]] + log(Rvec * Recruits[t])) + Nrand[ ,q] # add var. in recruitment
      Nf[Nf < 0] = 0
      
    
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
      avg.likel[t] = mean(ftmp[t,])
      avg.likel[which(!is.finite(avg.likel))] <- log(1e-320) # smallest non-infinite value
      
      ## Step 5: Resample    
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
  
    Pop_list[[y]] <- N
    Avg.likel[y] <- sum(avg.likel, na.rm = T) #save sum of each yr to take sum 
  } # end of year loop - Nact list
  
  
  Pfilter_data <- list(likelihood = -1 *sum(Avg.likel, na.rm = T),
                       Fit = Pop_list,
                       True = nact,
                       recruits = Recruits)
  
  return(Pfilter_data$likelihood)
}
# scale = 0.00001
# 