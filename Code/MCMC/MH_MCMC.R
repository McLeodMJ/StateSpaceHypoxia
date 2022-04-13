### MCMC - Metroplosis Hastings with Particle Filter likelihood 

#' Chain indicates which MCMC chain is currently being run (if not given, or given as NaN, 
#' then it will run the number of chains specified and run them in series (rather than parallel)

#setwd("~/Documents/THESIS/StateSpaceHypoxia/Code")
setwd("~/Box Sync/McLeod_thesis/Code")#Load MCMC files/ fxns
#load("../Data/DO_sims.Rdata")
#load("./Results/Priors_fromNO_Hypox.Rdata") # priors from MCMC bias runs [occ & det]?
source("./Library/logit.R")
source("./Library/inv.logit.R")
source("./Library/pois.likel.R")

# Load IPM Files/ fxns
load("./Results/Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors
#load("../Data/Hist_spp.Rdata") 
source("./Library/IPM.R")
source("./Library/params.R")
source("./Library/kernmat.R")
source("./Library/fecmat.R")
source("./Library/p.filter.R")
Sp_selex <- read.csv("../Data/Sp.selex.csv") 

library(matrixStats)
load("../Data/MH.MCMC.Rdata")

## Step 1: Load parameters
#det = 0.6
hypox_p = MH.dat$hypox_p
fi = MH.dat$fi
Rec_var = rep(MH.dat$Rec_var, MH.dat$time)
sigma_p = 1

# Step 2: prep MCMC

chains = 4
iters = 100000
# 100,000
# burn-ins manually 
# thinning - 1/ length(parm_vec)

mc_str = list(
  Post.prob = matrix(NA, ncol = ((iters - burn) * thin), nrow=chains),
  parm_vec = vector("list", chains),
  acc = rep(NA, chains),
  rej = rep(NA, chains),
  rej2 = rep(NA, chains)) # stores the rejection values

for(c in 1:chains){

   # initial Priors vector
   prior_mu = c(abs(MH.dat$prior_mu[3:4]), rep(MH.dat$prior_mu[5], MH.dat$time), MH.dat$prior_mu[6])  # hypox_p, fi, rec.variation, sigma_p (process) # inverse gamma hyperparameter for process error. Low confidence.
 # prior_mu = log(prior_mu) # CANNOT RUN B/C MOSTLY -INF
   parm_vec = matrix(NA, ncol= length(prior_mu), nrow=iters)
   parm_vec[1,] <- prior_mu
   prior_SD = c(MH.dat$prior_sd[3:4], rep(MH.dat$prior_sd[5], MH.dat$time), MH.dat$prior_sd[6])
   # prior_SD = log(prior_SD) # CANNOT RUN B/C WANT TO ALLIGN WITH MU
   CV = 1 / c(0.5, 1, 5:4, 3:40)  #may need to be a matrix
    stepmax = length(CV)
   
   # stores check-in values of MCMC 
   k = 1 # parameter counter
   m = 1 # iterations - delayed rejection counter 
   rej = 0 # count number of rejections
   rej2= 0  # count number of stepped rejections
   acc = 0 # count number of acceptances
   
   ## step 3: likelihood estimate
   # Run particle filter
  likel <- p.filter(MH.dat$Ninit, MH.dat$Nact, MH.dat$hypox, F, "Dsole", hypox_p, fi, MH.dat$mesh, MH.dat$Q, MH.dat$time, Rec_var, MH.dat$cv_q, sigma_p, 1)$likelihood
 
  
  # Prior calcs
  PriorL <- rep(NA, length(prior_mu))
  Post.prob <- NULL
  #PriorL[1] <- dnorm(det, prior_mu[1], prior_SD[1]) #det # remember to exp. afterwards 
  PriorL[1] <- dnorm(hypox_p, prior_mu[1],prior_SD[1]) #hypox_p
  PriorL[2] <- dnorm(fi, prior_mu[2],prior_SD[2]) #fi
  PriorL[3: c(length(PriorL) -1)] <- dnorm(Rec_var, prior_mu[3], prior_SD[3]) # rec_var
  PriorL[length(PriorL)] = sum(log(max(1e-26, dgamma(x=1 /parm_vec[1, ncol(parm_vec)] , shape = sigma_p, scale = (1 /prior_mu[length(prior_mu)]) /  (sigma_p-1)) ))) # process error
   # check so that never get value below zero
  Post.prob[1] <- sum(likel[2:length(likel)]) + sum(PriorL)

  # MCMC iterations
for(m in 2:iters){
  #m=1
    if((m%%1) == 0){ 
       #  matlab: if mod(m,1000) ==0
      print(paste('m =', as.character(m))) # counter - do we want to print the number for m as a character?
  }
  
      next_step = FALSE
      step = 1
      
      while(next_step == FALSE){ 
      
      # simulate candidate values
        # conjugates for random draws 
    cand_vec = parm_vec[m-1,]
   
      if (!is.na(prior_SD[k])){
        cand_vec[k] <- rnorm(1, cand_vec[k], prior_SD[k] * CV[step]) #det, hypox.p, fishing rate and recruitment (normal)
      }else{
        cand_vec[k] <- 1/rgamma(n=1, shape = sigma_p / CV[step], scale=c(1 / cand_vec[k] / (sigma_p / CV[step])) ) # process error - inverse gamma 
      }
      cand_vec[1:c(length(cand_vec) - 1)] = abs(cand_vec[1:c(length(cand_vec) - 1)]) # can't have a negative fishing rate or recruitment!
    
  P_filter_out <- p.filter(MH.dat$Ninit, MH.dat$Nact, MH.dat$hypox, F, "Dsole", cand_vec[1], exp(cand_vec[2]), MH.dat$mesh, MH.dat$Q, MH.dat$time, exp(cand_vec[3:c(length(cand_vec) - 1)]), MH.dat$cv_q, cand_vec[length(cand_vec)], 1)
 ## DONT THINK HYPOX.P HAS TO BE EXP.     
likel <- P_filter_out$likelihood
 
 # Prior calcs
       PriorL <- rep(NA, length(prior_mu))
     #PriorL[1] <- dnorm(cand_vec[1], prior_mu[1], prior_SD[1], log = TRUE) #det
     PriorL[1] <- dnorm(cand_vec[1],prior_mu[1],prior_SD[1], log = TRUE) # hypox
     PriorL[2] <- dnorm(cand_vec[2], prior_mu[2],prior_SD[2], log = TRUE) # fi
     PriorL[3:c(length(PriorL) - 1)] <- dnorm(cand_vec[3: c(length(PriorL) - 1)], prior_mu[3], prior_SD[3], log=TRUE) # recr
    
     PriorL[length(PriorL)] = log(max(1e-26, dgamma(x=1 /cand_vec[ncol(parm_vec)] , shape = sigma_p, scale = (1 /prior_mu[length(prior_mu)]) /  (sigma_p-1)) )) # process error
     L_cand <- sum(likel[2:length(likel)]) + sum(PriorL) # Candidate likelihood
  
  ## Metropolis-Hastings step
  if(!is.na(L_cand) & !is.infinite(L_cand)){
    MH_prob <- pmin(1, exp(L_cand - Post.prob[m-1])) # NEED TO MAKE SURE THE POST.PROB & L_CAND ARE ALLIGNED
    # IN ORDER FOR THIS TO RUN PROPERLY. CURRENTLY GETTING: INF -> Mh_prob = 1 -> K = true
  }else{
    MH_prob = 0
  }
  # check to see if this value is larger than a random value btw 0-1
      Acc <- MH_prob > runif(1)
 
  # Delayed rejection step 
      if(Acc == TRUE){ # accept proposal 
      parm_vec[m,] = cand_vec
      Post.prob[m] = L_cand
      next_step = TRUE
      step = 1
      acc = acc + 1
      
      } else if(Acc == F & step < stepmax){ # move to alternate proposal
        step = step + 1
        rej2 = rej2 + 1
      
      } else{ # if actually reject
        parm_vec[m,] = parm_vec[m-1,]
        Post.prob[m] = Post.prob[m-1]
        next_step = TRUE
        rej = rej + 1
        } # end if K
      
      } # end if delayed rejection - while loop
        
    # advance counters
      k = k+1
      if(k > length(cand_vec)){
        k = 1 } # after the while loop is broken we need to move onto next param 
   } # end loop over M iterations
  
  ### Burn in and thinning
  Post.prob <- Post.prob[c(burn + 1) : length(Post.prob) ] * thin

  mc_str$Post.prob[c,] = Post.prob
  mc_str$parm_vec[[c]] = parm_vec
  mc_str$acc[c] = acc
  mc_str$rej[c] = rej
  mc_str$rej2[c] = rej2
    
} # end of chains loop 
  
#save(mc_str, file= "./Results/MH.MCMC_Results.Rdata")  
