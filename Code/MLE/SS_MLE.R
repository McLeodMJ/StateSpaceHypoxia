
source("logit.R")
source("inv.logit.R")
source("pois.likel.R")

# Load IPM Files/ fxns
source("IPM.R")
source("params.R")
source("kernmat.R")
source("fecmat.R")
source("p.filter.R")
#source("MH.MCMC.R")
Sp_selex <- read.csv("Sp.selex.csv") 
load("Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors
load("Dsole_true_params_log.Rdata") # log of the typical Dsole data
#load("Nact_DO.Rdata")

# library(matrixStats)
# library(parallel)
 library(foreach)
library(doParallel)
 library(invgamma)
library(bbmle)
library(tidyverse)
library(dplyr)
#library(mcreplicate)
library(httr)
library(broom)

pars <- params("Dsole")

# test parameters over MCMC
test_params <- cbind( c(MH.Sim.Dsole$hypox_p, rep(MH.Sim.Dsole$hypox_p, 2), 0.5, 8) ,#hypox_p
                      c(log(pars$f), log(0.001), log(0.2), rep(log(pars$f), 2) )) # fi


#params_dlog <- data.frame( Fi = dnorm(test_params[,1], MH.Sim.Dsole$prior_mu[2], MH.Sim.Dsole$prior_sd[2], log=T), # fi
#                          Rec_vec = dnorm(0.1, MH.Sim.Dsole$prior_mu[3], MH.Sim.Dsole$prior_sd[3], log=T),
#                        Sigmpa_p = dnorm(0.1,  MH.Sim.Dsole$prior_mu[4], MH.Sim.Dsole$prior_mu[4], log = T))
#
## parameters to test out on lognorm scale using rlnorm()
params_log <- data.frame( Fi = c(log(MH.Sim.Dsole$prior_mu[2]), log(0.0001), log(0.2), log(MH.Sim.Dsole$prior_mu[2]), log(MH.Sim.Dsole$prior_mu[2])), # fi
                           hypox_p = c(4, 4, 4, 1, 8),
                          sigma_p = c(log(0.001), log(0.01), log(0.1), log(0.5), log(2)))

i=1
log_start <- list(hypox_p = params_log[1,2],
                         fi = params_log[1,1],
                         cv_q = 0.1,
                         sigma_p = params_log[i,3],
                         rec1 =  log(0.1),
                         rec2 =  log(0.1),
                         rec3 =  log(0.1),
                         rec4 =  log(0.1),
                         rec5 =  log(0.1),
                         rec6 =  log(0.1),
                         rec7 =  log(0.1),
                         rec8 =  log(0.1),
                         rec9 =  log(0.1),
                         rec10 = log(0.1),
                         rec11 = log(0.1),
                         rec12 = log(0.1),
                         rec13 = log(0.1),
                         rec14 = log(0.1),
                         rec15 = log(0.1),
                         rec16 = log(0.1),
                         rec17 = log(0.1),
                         rec18 = log(0.1),
                         rec19 = log(0.1),
                         rec20 = log(0.1))


low <-c("hypox_p" = -100, "fi" = -100, "sigma_p"= -100, "rec1"= -100, "rec2"= -100, "rec3"= -100, "rec4"= -100, "rec5"= -100, "rec6"= -100, "rec7"= -100, "rec8"= -100, "rec9"= -100, "rec10"= -100, "rec11"= -100,"rec12"= -100, "rec13"= -100, "rec14"= -100, "rec15"= -100, "rec16"= -100, "rec17"= -100, "rec18"= -100, "rec19"= -100, "rec20"= -100)
upp <-c("hypox_p" = 100, "fi" = 100, "sigma_p"= 100, "rec1"= 100, "rec2"= 100, "rec3"= 100, "rec4"= 100, "rec5"= 100, "rec6"= 100, "rec7"= 100, "rec8"= 100, "rec9"= 100, "rec10"= 100, "rec11"= 100,"rec12"= 100, "rec13"= 100, "rec14"= 100, "rec15"= 100, "rec16"= 100, "rec17"= 100, "rec18"= 100, "rec19"= 100, "rec20"= 100)

load("./Results/SS_fi_std.Rdata")
ests <- coef(SS_fi_std)
ests <- ests[c(1:2, 4:length(ests))]

# neg. LL for fxn that contains pre & post-molt LL's
start_time <- Sys.time()
print(start_time)

SS_fi_std <- mle2(minuslogl =  p.filter, #p.filter.wcbgts
                 start = log_start,
                 fixed = list(cv_q = 0.1),
                 data = list(dat = MH.Sim.Dsole.log,
                             scale = 1),
                 method = "L-BFGS-B",
                  lower = low, 
                 upper = upp,
                 control=list(maxit=5000, trace=TRUE))
                 
               # skip.hessian = T,
               # control = list(maxit = 50),
               # trace = T)
    
save(SS_log_final5, file ="./Results/SS_log_final5.Rdata")
end_time <- Sys.time()
end_time - start_time


test <- mle2(minuslogl =  p.filter, #p.filter.wcbgts
                        start = log_start,
                        fixed = list(cv_q = 0.1),
                        data = list(dat = MH.Sim.Dsole.log,
                                    scale = 1))
save(test, file="test.Rdata")


# method = "L-BFGS-B", 
# lower = rep(-Inf, 23), upper = rep(Inf, 23),
# control = list(trace = 5, fnscale = -1) ) 


# control = list(1000) $ max iterations
               #  method = "L-BFGS-B", 
                # lower = rep(-Inf, 23), upper = rep(Inf, 23),
                # control = list(trace = 5, fnscale = -1) ) )

#MLE_Results_Total[[t]] <- summary(state.space.mll)
#}
 
 
 #mclapply(p.filter, mle2, mc.cores = 2, start = start_list, fixed = fixed_list, data = data_list, method = "L-BFGS-B",  lower = rep(-Inf, 23), upper = rep(Inf, 23), control = list(trace = 5, fnscale = -1))
# summary(state.space.mll) 
# SS_results <- summary(state.space.mll) 
 #save(MLE_Results_Total, file= "./Results/MLE_Results.Rdata")
 
 setwd("~/Documents/THESIS/StateSpaceHypoxia/Code")
 #setwd("~/Box Sync/McLeod_thesis/Code")#Load MCMC files/ fxns
 load("../Data/DO_sims.Rdata")
# load("./Results/Priors_fromNO_Hypox.Rdata") # priors from MCMC bias runs [occ & det]?
 source("./Library/logit.R")
 source("./Library/inv.logit.R")
 source("./Library/pois.likel.R")
 
 # Load IPM Files/ fxns
 load("./Results/Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors
 #load("../Data/Hist_spp.Rdata") 
 #load("../Data/MH.Sim.Dsole.Rdata")
# load("./Library/Dsole_true_params_log.Rdata")
 source("./Library/IPM.R")
 source("./Library/params.R")
 source("./Library/kernmat.R")
 source("./Library/fecmat.R")
 source("./Library/p.filter.R")
 source("./Library/p.filter.WCGBTS.R")
 #source("./Library/MH.MCMC.R")
 Sp_selex <- read.csv("../Data/Sp.selex.csv") # 

 
 # library(matrixStats)
 # library(parallel)
 library(foreach)
 library(doParallel)
 library(invgamma)
 library(bbmle)
 library(tidyverse)
 library(dplyr)
 #library(mcreplicate)
 library(httr)
 library(broom)
 
 pars <- params("Dsole")
 
 # test parameters over MCMC
 test_params <- cbind( c(MH.Sim.Dsole$hypox_p, rep(MH.Sim.Dsole$hypox_p, 2), 0.5, 8) ,#hypox_p
                       c(log(pars$f), log(0.001), log(0.2), rep(log(pars$f), 2) )) # fi
 
 
 i=1
 log_start <- list(hypox_p = test_params[i,1],
                   fi = test_params[i,2],
                   cv_q = 0.1,
                   sigma_p = log(0.5),
                   rec1 =  log(0.1),
                   rec2 =  log(0.1),
                   rec3 =  log(0.1),
                   rec4 =  log(0.1),
                   rec5 =  log(0.1),
                   rec6 =  log(0.1),
                   rec7 =  log(0.1),
                   rec8 =  log(0.1),
                   rec9 =  log(0.1),
                   rec10 = log(0.1),
                   rec11 = log(0.1),
                   rec12 = log(0.1),
                   rec13 = log(0.1),
                   rec14 = log(0.1),
                   rec15 = log(0.1),
                   rec16 = log(0.1),
                   rec17 = log(0.1),
                   rec18 = log(0.1),
                   rec19 = log(0.1),
                   rec20 = log(0.1))
 
 low <-c("hypox_p" = -Inf, "fi" = -Inf, "sigma_p"= -Inf, "rec1"= -Inf, "rec2"= -Inf, "rec3"= -Inf, "rec4"= -Inf, "rec5"= -Inf, "rec6"= -Inf, "rec7"= -Inf, "rec8"= -Inf, "rec9"= -Inf, "rec10"= -Inf, "rec11"= -Inf,"rec12"= -Inf, "rec13"= -Inf, "rec14"= -Inf, "rec15"= -Inf, "rec16"= -Inf, "rec17"= -Inf, "rec18"= -Inf, "rec19"= -Inf, "rec20"= -Inf)
 upp <-c("hypox_p" = Inf, "fi" = Inf, "sigma_p"= Inf, "rec1"= Inf, "rec2"= Inf, "rec3"= Inf, "rec4"= Inf, "rec5"= Inf, "rec6"= Inf, "rec7"= Inf, "rec8"= Inf, "rec9"= Inf, "rec10"= Inf, "rec11"= Inf,"rec12"= Inf, "rec13"= Inf, "rec14"= Inf, "rec15"= Inf, "rec16"= Inf, "rec17"= Inf, "rec18"= Inf, "rec19"= Inf, "rec20"= Inf)
 
 
 
 start_time <- Sys.time()
 print(start_time)
 
 Final_SS_std <- mle2(minuslogl =  p.filter, #p.filter.wcbgts
                         start = log_start,
                         fixed = list(cv_q = 0.1),
                         data = list(dat = MH.Sim.Dsole.log_standard,
                                     scale = 1),
                         method = "L-BFGS-B",
                         lower = low, 
                         upper = upp,
                         control=list(maxit=5000, trace=TRUE))
 
 save(Final_SS_std, file ="./Results/Final_SS_std.Rdata")
 end_time <- Sys.time()
 end_time - start_time
 
