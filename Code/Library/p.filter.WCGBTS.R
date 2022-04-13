####################################
# particle filter function
#' dat: simulated IPM data --> equilibrium
#' Nact: actual data from the WCGBTS [if checking if p.filter works - then this is 2nd run of IPM]
#' selc: F means we are running simulations and do NOT need the WCGBTS selc- T means we need WCGBTS selc. 
#' fish: "Dsole", "Lingcod", "Yrock", "Grock" to call param for species of interest
#' fi: fishing rate - take fishing rate found from SPR Analysis OR pars$f generally fits to call from param() fxn
#' hypox_a: dep. hypoxia parameter - based on 1.43 threshold for hypoxia - hypoxa=7/1.43 --> hypox_a= 4.9
#' hypox_b: intercept hypoxia parameter - if DO = 0 & prob is very low [0.001] logit(0.001/.999) --> hypox_b = -7
#' mesh: mesh size - defines nrows in matrix
#' Q: no. of particles 
#' time: running iterations to this time
#' rec_var: variation for recruitment variation
#' cv_q: random variation introduced for each particle 
#' sigma_p: variation for process error
#' scale: used to reduce pop sizes down to more realistic values per trawl - best scale = 1e-6 multiplied by Nf

#######################################

p.filter.WCGBTS <- function(dat, hypox_a, hypox_b, fi, cv_q, sigma_p, scale, rec1, rec2, rec3, rec4, rec5, rec6){
  # return data from list
  nact = dat$N.act
  hypox = dat$DO
  fish = dat$fish
  mesh = dat$Mesh
  Q = dat$Q
  rec_var <- list("2010" = rec1, 
                     "2011" = rec2,
                     "2012" = rec3, 
                     "2013" = rec4, 
                     "2014" = rec5,
                     "2015" = rec6)
  
# time must be formatted by matrices of each year (ncol - trawl sites)
   time = length(nact)+1 
  
### Step 1:  create Inital SAD using IPM && fi parameter being estimated by p.filter
    N0 <- IPM(fish, fi, 100, mesh, NA, rep(0, 100))$SAD
    # format data to be scaled down to more realistic sizes
  
    ############### Option 1 ###################   
    # N0 <- N0 * scale #(make starting values smaller)
    # BUT then do we scale down the R0 and S0 values used in Beverton-Holt? (line - 126) anywhere else?
    ############### Option 2 ###################   
   # nact / scale # make actual values larger
    # BUT then do we need to make cv_q larger? (variation between particles for intial values) 
    # what about sigma_p? (for adding process error)
 
       #update progress 
    start_time <- paste('Start time: ', Sys.time(), sep = '')
    write(start_time, file = 'MLE_WCGBTS_progress.txt', append = TRUE)
    
 
## Step 2: Initialize the process model 
    N <- matrix(NA, nrow = mesh, ncol = time)
    
    # Generate Q particles (independent simulations of N):
    Nf <- array(rep(NA, mesh * time * Q), c(mesh, time, Q)) #pop.dist. x time x particles
    Nf[,1,] <- N0 #inputs the last column of pop. equil. data
    
    # Simulate random variation for first particle
    part_sims <- rnorm(length(Nf[,1,]), 0, cv_q)
    dim(part_sims) <- dim(Nf[,1,]) 
    Nf[,1,] <- Nf[,1,] + part_sims
    Nf[Nf<0] = 0 # check for all values greater than zero 
    
    
## Step 3: Initialize the Operating model 
  # resample for accuracy
    avg.likel <- rep(NA, (time-1)) #average of the ftmp
    #likel <- matrix(NA, nrow=mesh, ncol= Q)
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
    pars <- params(fish) # calls parameters from SS
  
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
    
    # Factor in WCGBTS size-Selectivity
     WClen <- pars$selc 
     
## Step 5: run model through time T i.e. years
    Recruits <- NULL
    
   for (t in 2:time){
     # Detection parameter dependent on hypoxia parameter and DO data
     hyp <- hypox[[t-1]]$mean_hyp
     det <- inv_logit(hyp * hypox_a + hypox_b)
     det[det < 1e-6] = 1e-6
     det[det > (1-1e-6)] = (1-1e-6)
     
     
     #population by year
    Nact <- nact[[t-1]]
   
     pois_ll <- array(rep(0, mesh*ncol(Nact)*Q), c(mesh, ncol(Nact), Q)) # stores pois loglik across obsv /time
     Ftmp <- matrix(NA, ncol(Nact), Q)
     Ftmp[1,] <-rep(1,Q)
    
 ## Step 6: Advance the model, one particle at a time
    for(q in 1:Q){
      # integration
      Nf[,t,q] <- K %*% Nf[,t-1,q]  * dx  # midpoint rule integration
      E <- sum(Fe * Nf[,t-1, q] * dx )
      
      # Beverton-holt density-depedence
      Recruits[t] <- as.vector((4 * pars$steep * exp(pars$R0) * E)/ ((pars$S0 * (1 - pars$steep)) + (E * (5 * pars$steep - 1))) ) 
      # THIS^^ gives a reasonable value but does it make sense to apply this way??? - as.vector((4 * pars$steep * exp(pars$R0*scale)* E)/ ((pars$S0*scale * (1 - pars$steep)) + (E * (5 * pars$steep - 1))) ) 
      Nf[,t, q] = Nf[,t, q] +  exp(rec_var[[t-1]] + log(Rvec * Recruits[t])) + Nrand[ ,q] # add var. in recruitment
      Nf[Nf < 0] = 0
      
       
## Step 7: pois-likelihood        
  ############################### LIKELIHOOD ############################
          
    ## additional for loop to run over each trawl sites 
    for(oo in 1:ncol(Nact)){
        pois_ll[,oo,q] <- dpois(x=Nact[,oo], lambda = (Nf[,t,q]/ exp(scale)) * WClen * dx , log=T)
           #pois_mle <- data.frame(lambda_vals = Nf[,t,q] * dx * WClen * scale) %>%
           #  rowwise() %>%
           #   mutate(log_likelihood = pois.likel(y = Nact[,oo], lambda = lambda_vals)) %>%
           #   ungroup() # avoids issues with dpois() over vectorized data
           # pois_ll[,oo,q] <- pois_mle$log_likelihood
       } # end of oo loop
     } # end of Q loop

      like.min <- min(pois_ll[is.finite(pois_ll)], na.rm = T) 
      #add ifelse statement for when bounds are reached and entire likel matrix returns NaNs, like.min will return Inf
      if(like.min == Inf){
        like.min = -1e50
      }
      pois_ll[which(!is.finite(pois_ll))] <-  like.min 
      
## Step 8: zero-inflated model with detectability as a fxn of hypoxia 
    for(o in 1:ncol(Nact)){
        ## factoring in detection parameter to likelihood data
        if(sum(Nact[ ,o]) == 0){ 
          # absent
         # Ftmp[t,] <- colSums(likel,na.rm=TRUE) * c(1 - det[t]) / c(sum(det) * time)
          Ftmp[o,] <- colSums(pois_ll[,o,],na.rm=TRUE) * c(1 - det[o]) / c(sum(det) * ncol(Nact))
        }else{
          #present
         # Ftmp[t,] <- colSums(likel,na.rm=TRUE) * det[t] /c(sum(det) * time)
          Ftmp[o,] <- colSums(pois_ll[,o,],na.rm=TRUE) * det[o] /c(sum(det) * ncol(Nact))
        }
    } # end of observation loop
        ftmp[t,] <- colSums(Ftmp) #sum over particles for all observations
   avg.likel[t-1] = mean(ftmp[t,])
     
## Step 9: Resample / Particle filter: (following Knape & deValpine 2012)  
      #Wgt = cumsum(Ftmp[o,]/ sum(Ftmp[o,]))
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
    
    ## sum of the negative log likelihood averaged across years
    NLL <- -1 *sum(avg.likel, na.rm = T)
    
    # get updates of the parameters in real time 
    update <- paste('// hypox_a=',hypox_a, '/ hypox_b=', hypox_b, '/ fi=', fi, '/ sigma=', sigma_p, '/ like.min=', like.min,'/ NLL=', NLL, sep = ' ')
    write(update, file = 'MLE_WCGBTS_progress.txt', append = TRUE)
    
  Pfilter_data <- list(likelihood = NLL,
                       Fit = N,
                       True = nact,
                       recruits = Recruits)
  
  return(Pfilter_data$likelihood)
}
# scale = 0.00001
 