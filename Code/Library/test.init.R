########################################################################################
### Function to run multiple MCMCs to test different intial values for parameters
# dat = good_yr or bad_yr (DO simulations - drawn from 'Data' folder)
#########################################################################################

test.init <- function(dat, priors, tru) {
  start_time <- paste('Start time: ', Sys.time(), ' - ', dat, ' - ', tru, sep = '')
  write(start_time, file = 'progress.txt', append = TRUE)
  
  N = ncol(dat) #N=length(dat) if solo try
  simz_list <- list()
  model <- stan_model(occ.mod1)
  
  #default return as list
  #foreach (j = 1:nrow(tru)) %dopar% {
  for (j in 1:nrow(tru)){
    
    occ = tru$occ[j]
    det = tru$det[j]
    hypox.p = tru$hypox.p[j]
    
    # df for timeseries simulations 
    simz <- data.frame(occ = rep(occ, nrow(dat)),
                       det = rep(det, nrow(dat)),
                       hypox.p = rep(hypox.p, nrow(dat)),
                       occ_post = rep(0, nrow(dat)), 
                       det_post = rep(0, nrow(dat)),
                       hypox.p_post = rep(0, nrow(dat)),
                       occ_bias = rep(0, nrow(dat)),
                       occ_sd = rep(0, nrow(dat)),
                       det_bias = rep(0, nrow(dat)),
                       det_sd = rep(0, nrow(dat)),
                       hypox.p_bias = rep(0, nrow(dat)),
                       hypox.p_sd = rep(0, nrow(dat)) )
    #solo simz df
   #simz <- data.frame(occ = rep(occ,1),det = rep(det,1),hypox.p = rep(hypox.p, 1),occ_post = rep(0, 1), det_post = rep(0, 1),hypox.p_post = rep(0,1),occ_bias = rep(0, 1)),occ_sd = rep(0,1),det_bias = rep(0, 1),det_sd = rep(0,1),hypox.p_bias = rep(0, 1),hypox.p_sd = rep(0, 1) )
    #################################### factoring in hypoxia as a covariate ############################################
    for (s in 1:nrow(dat)){
      occ.full <- inv_logit(logit(occ) + hypox.p * dat[s,]) #solo = dat
      
      # encounters [0 or 1]
      enc <- rbinom(N, 1, occ.full) * rbinom(N, 1, det)
      
 ###############################  Creating data in list format necessary for Stan ####################################
      occ.stan.dat <- list(N = N,
                           enc = enc,
                           hypox = dat[s,],  #solo = dat
                           prior_mu = c(priors$mean[1], priors$mean[2], 0),
                           prior_sd = c(priors$sd[1], priors$sd[2], 10)) #flatter- 1
      
      
  ###############################################  STAN MODEL  ############################################ 
#occ.stan <- stan(model_code = Occ.mod,
#           pars = c('occ','det','hypox_p'),
#           data = occ.stan.dat,
#           chains = 4,
#           cores = 4,
#           warmup = 2000,
#           iter = 8000,
#           control = list(adapt_delta = 0.95))
     
      fit <- sampling(object = model, 
                      data = occ.stan.dat,
                      chains = 4, 
                      cores = 4,
                      iter = 8000, 
                      warmup = 2000,
                      control = list(adapt_delta = 0.99))
   
      
      sum <- as.data.frame( summary(fit)$summary )
      
      # my janky for loop to check if rhat is too large or effect size too small 
      x <- rep(0, nrow(sum))
      for(n in 1: nrow(sum)){
        # insufficient rhat >=1.05 & effect size <100
        if(sum$Rhat[n] >= 1.05 | sum$n_eff[n] <= 100) {
          x[n] <- NA
          } }
        x <- na.omit(x) # adjusts length for next ifelse statement 
        
       if(length(x) < nrow(sum)) {
          #insufficient 
          simz[s,4:6] <- rep(NA, 3)
        }else{
          #sufficient 
          simz[s,4:6] <- c(sum$mean[4], sum$mean[5], sum$mean[3])
        }
      
      # calculate bias
      simz$occ_bias[s] <-  (simz$occ_post[s] - occ) / simz$occ_post[s]
      simz$det_bias[s] <- (simz$det_post[s] - det) / simz$det_post[s]
      simz$hypox.p_bias[s] <- (simz$hypox.p_post[s] - hypox.p) / simz$hypox.p_post[s]
      
      # retrieve the SD 
      simz$occ_sd[s] <- sum$sd[4]
      simz$det_sd[s] <- sum$sd[5]
      simz$hypox.p_sd[s] <- sum$sd[3]
      
    } #end s for loop
    simz_list[[j]] <- simz 
    
    
  } #end j for loop
  if (j %% 2 == 0) {
    update <- paste(Sys.time(), ' - ', dat, ' - ', j, 'done!', sep = '')
    write(update, file = 'progress.txt', append = TRUE)
  }
  
  save(simz_list, file= "Simz_list.Rdata")
  #save(simz_list, file= paste('~/Documents/MCMC/Library/Simz_listzzzzz', name, '.Rda', sep=''))

} #end fxn

