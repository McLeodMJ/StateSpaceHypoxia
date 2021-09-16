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

# load("../Data/Hist_spp.Rdata")
# load("./Results/Expected_Fishing_rate.Rdata")
# source("./Library/params.R")
# source("./Library/kernmat.R")
# source("./Library/fecmat.R")
# Sp_selex <- read.csv("../Data/Sp.selex.csv") 


#######################################

p.filter <- function(dat, Nact, hypox, fish, hypox_p, fi, mesh, Q, time, rec_var, cv_q, sigma_p){
  
  # Particle filter: (following Knape & deValpine 2012)
  N <- matrix(NA, nrow = nrow(dat), ncol = time)
  
  # Generate Q particles (independent simulations of N):
  Nf <- array(rep(NA, nrow(N) * time * Q), c(nrow(N), time, Q)) #pop.dist. x time x particles
  Nf[,1,] <- dat[,time] #saves last column of equil. data
  
  # Simulate random variation for first particle
  Nf[,1,] <- Nf[,1,] + rnorm(nrow(N), 0, cv_q * dat[,time])
  Nf[Nf<0] = 0 # check for all values greater than zero 
  
  # Step 1: Initialize the model 
    # resample for accuracy
    ftmp <- matrix(NA, nrow= time, ncol = Q)
    Avg.likel <- rep(0, time) #average of the ftmp
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
  
  ## Step 2: simulate encounters for Null or absences based on hypoxia parameter 
  
  #threshold -50
  #occ.full.2 <- det + inv_logit(hypox * hypox_p) #linear dep. on hypoxia
  # may want to replace data with simulated pop and set a threshold
  # threshold is blanket greater .5 then it is true. [check the hypox for if larger is 0 or 1]
  
  # for loop but of occ.full.2 
  #for(i in 1:time){
   # if(occ.full.2[i] <= 0.89){
    #  ftmp[i,] <- rep(NA, Q) # think its due to hypox replace zero observ. moments with -999 [true zero =]
    #}}
  
  ## Step 3: run model through time T 
    # Create the kernel:
    pars <- params(fish)
    
    meshmax = pars$Linf * 2
    x <- seq(from = 1, to = meshmax, length.out = mesh)
    dx = diff(x)[1]
    
  ## size distribution of recruits
  Rvec <- dnorm(x, pars$Rlen, pars$Rlen.sd)  
  
  K <- kernmat(x, pars, fi)
  Fe <- fecmat(x, pars)
  
  for (t in 2:time){
    
    # Advance the model, one particle at a time
    for(q in 1:Q){
      Nrand = matrix(rnorm(mesh * Q, 0 , sigma_p), mesh, Q) #processs error
      
      # integration
      Nf[,t,q] <- K %*% Nf[,t-1,q]  * dx  # midpoint rule integration
      E <- sum(Fe * Nf[,t-1,q]) * dx 
      
      # Beverton-holt density-depedence
      Recruits <- as.vector((4 * pars$steep * exp(pars$R0) * E) / ((pars$S0 * (1 - pars$steep)) + (E * (5 * pars$steep - 1))))
      # Add variation in recruitment 
      
      #RR <- exp(rnorm(1, mean = log(Recruits) - ((cv * log(Recruits))^2)/2, sd= cv * log(Recruits) ))# change cv to 0 for no variation
      
      Nf[,t, q] = Nf[,t, q] + rec_var[t] * (Rvec * Recruits) + Nrand[ ,q]  # move the model forward for each particle 
    #  Nf(:,t,q) = kmat*Nf(:,t-1,q) + RR + Nrand(:,:,q);  % Eq 2 from White et al. PLoS ONE
    
      } # end of Q loop
    Nf[Nf < 0] = 0
    
    # WCGBTS Selectivity
      WClen <- pars$selc
      
      ############################### LIKELIHOOD ############################
      det <- inv_logit(hypox * hypox_p)
      
      likel <-dpois(Nact[,t], Nf[,t,] * dx * WClen, log= T)
      #= sum(log(max(realmin,poisspdf(repmat(Nact(:,t),[1,1,Q]),(NT(t-length(Tpre))*dy*Nf(:,t,:)).*repmat(OKlen(:),[1,1,Q])))));
      likel[which(!is.finite(likel))] <- 0
      # eqn we worked through on 9/1
      #ftmp[t,] = colSums(likel) * occ.full.2 + (1 - occ.full.2) * (sum(Nact[,t]) == 0)  #log-likelihood vector
      ### eqn 1
      #ftmp[t,] = ifelse(sum(Nact[,t]) == 0, (1 - occ.full.2), colSums(likel) * occ.full.2 )  #log-likelihood vector
      ### eqn 2/3
      #likel_b <- dbinom(x = Nact[,t], size = Nf[,t,], prob = det[t])
       if(sum(Nact[,t]) == 0){ 
         # absent
        ftmp[t,] <- colSums(likel) * c(1 - det[t])
        }else{
          #present
         ftmp[t,] <- colSums(likel) * det[t]
        }
      Avg.likel[t] = mean(ftmp[t,])
      
      ### alternative LLs
      #likel_b = c(log( det[t] + exp(log(1 - det[t]) + 1)))
      #likel_p = log(1 - det[t]) + colSums(dpois(Nact[,t], Nf[,t,] * dx * WClen, log= T))
     
    
      # check - put in true parameters values from process model and check that it gives a higher likelihood for those over the recruitment was half of that year, fi- is 50% higher 
      
      
    # if negative binomial would be better THAN POISON
      #likel = dnbinom(x = Nact[,t], size = occ.full.2, mu = Nf[,t,] * dx * WClen, log= T)
      #size is target for no. of successful trials ???
      # cannot do p instead of size 
      
 
   ## Step 4: Resample    
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
  

  Pfilter_data <- list(likelihood = Avg.likel,
                         Pop = N)
  
  return(Pfilter_data)
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

 