
source("logit.R")
source("inv.logit.R")
source("pois.likel.R")

# Load IPM Files/ fxns
source("IPM.R")
source("params.R")
source("kernmat.R")
source("fecmat.R")
source("p.filter.R")
source("p.filter.WCGBTS.R")
source("obsv.sims.R")
#source("MH.MCMC.R")
Sp_selex <- read.csv("Sp.selex.csv") 
load("Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors
load("Dsole_true_params_log.Rdata") # log of the typical Dsole data
load("Nact_DO.Rdata")
load("Nact_site_yr.Rdata")

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
log_start <- list(hypox_p = test_params[1,1],
                  fi = test_params[1,2],
                  cv_q = 0.1,
                  sigma_p = log(0.5),
                  rec1 =  log(0.1),
                  rec2 =  log(0.1),
                  rec3 =  log(0.1),
                  rec4 =  log(0.1),
                  rec5 =  log(0.1))
                 

#low <-c("hypox_p" = -Inf, "fi" = -Inf, "sigma_p"= -Inf, "rec1"= -Inf, "rec2"= -Inf, "rec3"= -Inf, "rec4"= -Inf, "rec5"= -Inf, "rec6"= -Inf, "rec7"= -Inf, "rec8"= -Inf, "rec9"= -Inf, "rec10"= -Inf, "rec11"= -Inf,"rec12"= -Inf, "rec13"= -Inf, "rec14"= -Inf, "rec15"= -Inf, "rec16"= -Inf, "rec17"= -Inf, "rec18"= -Inf, "rec19"= -Inf, "rec20"= -Inf)
#upp <-c("hypox_p" = Inf, "fi" = Inf, "sigma_p"= Inf, "rec1"= Inf, "rec2"= Inf, "rec3"= Inf, "rec4"= Inf, "rec5"= Inf, "rec6"= Inf, "rec7"= Inf, "rec8"= Inf, "rec9"= Inf, "rec10"= Inf, "rec11"= Inf,"rec12"= Inf, "rec13"= Inf, "rec14"= Inf, "rec15"= Inf, "rec16"= Inf, "rec17"= Inf, "rec18"= Inf, "rec19"= Inf, "rec20"= Inf)


# neg. LL for fxn that contains pre & post-molt LL's
start_time <- Sys.time()
print(start_time)

test_DSole <- mle2(minuslogl =  p.filter.WCGBTS, #p.filter.wcbgts
                   start = log_start,
                   fixed = list(cv_q = 0.1),
                   data = list(dat = MH.Sim.Dsole.log_standard,
                               nact = Dsole.OR.WA,
                               hypox = Nact_DO,
                               scale = 1))
                   
                  #method = "L-BFGS-B",
                  #lower = low, 
                  #upper = upp,
                  #control=list(maxit=5000, trace=TRUE)

print(summary(test_DSole))
end_time <- Sys.time()
end_time - start_time

save(test_DSole, file="test_DSole.Rdata")





SS_fi_std <- mle2(minuslogl =  p.filter, #p.filter.wcbgts
                  start = log_start,
                  fixed = list(cv_q = 0.1),
                  data = list(dat = MH.Sim.Dsole.log,
                              scale = 1),
                  method = "L-BFGS-B",
                  lower = low, 
                  upper = upp,
                  control=list(maxit=5000, trace=TRUE))


save(SS_log_final5, file ="./Results/SS_log_final5.Rdata")

