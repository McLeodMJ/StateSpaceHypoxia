########################################################################################
### Function to run multiple MCMCs to test different intial values for parameters
# dat = good_yr or bad_yr (DO simulations - drawn from 'Data' folder)
#########################################################################################

test.init <- function(dat, tru) {
  start_time <- paste('Start time: ', Sys.time(), ' - ', dat, ' - ', 
                      tru, sep = '')
  write(start_time, file = 'progress.txt', append = TRUE)
  
  N = ncol(dat)
  simz_list <- list()
  
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
    
    #################################### factoring in hypoxia as a covariate ############################################
    for (s in 1:nrow(dat)){
      occ.full <- inv_logit(logit(occ) + hypox.p * dat[s,])
      
      # encounters [0 or 1]
      enc <- rbinom(N, 1, occ.full) * rbinom(N, 1, det)
      
      ###############################################  STAN MODEL  ############################################ 
      occ.mod2 <- "
    	data{
    		int<lower=1> N; //number of samples
    		int<lower=0> enc[N]; //encounters
    		vector[N] hypox; //hypoxia
    		vector[3] prior_mu;
    		vector[3] prior_sd;
    	}
    	parameters{
    		real logit_occ;
    		real logit_det;
    		real hypox_p; //hypoxia slope
    	}
    	transformed parameters{
    		real<lower=0,upper=1> occ; // occupancy
    		real<lower=0,upper=1> det; // detection
    		occ = inv_logit(logit_occ);
    		det = inv_logit(logit_det);
    	}
    	model{
    	  // storage
    	    vector[N] occ_eff;
    
    	  // weakly informative priors
    		logit_occ ~ normal(prior_mu[1],prior_sd[1]);
    		logit_det ~ normal(prior_mu[2],prior_sd[2]);
    		hypox_p ~ normal(prior_mu[3],prior_sd[3]);
    
    	// likelihoods
    		for(i in 1:N){
    			occ_eff[i] = logit_occ + hypox_p*hypox[i];
    			if(enc[i] > 0){
    				//the site was occupied, and you detected it
    				target += log_inv_logit(occ_eff[i]) + bernoulli_logit_lpmf(1| logit_det );
    				
    			}else{
    				//the site was occupied but you didn't detect it & the site was unoccupied
    				target += log_sum_exp(log_inv_logit(occ_eff[i]) + bernoulli_logit_lpmf(0| logit_det),       log1m_inv_logit(occ_eff[i]));
    			}
    		}
    	}"
      # target as a vector then take sum of all values for occupied and non-occupied. 
      
      ###############################  Creating data in list format necessary for Stan ####################################
      occ.stan.dat <- list(N = N,
                           enc = enc,
                           hypox = dat[s,],
                           prior_mu = c(0, logit(det), 0),
                           prior_sd = c(10, 10, 10)) #flatter- 1
      
      occ.stan <- stan(model_code = occ.mod2,
                       pars = c('occ','det','hypox_p'),
                       data = occ.stan.dat,
                       chains = 4,
                       warmup = 2000,
                       iter = 6000,
                       control = list(adapt_delta = 0.95)) 
      sum <- summary(occ.stan)
      
      
      # my fxn to check if rhat is too large or effect size too small 
      x <- rep(0, nrow(sum$summary))
      x <- convr.chk(sum$summary, x)
      
      # convergence detected -> report as NA   
      if(is.na(x)){
        simz[s,4:6] <- rep(NA, 3)
      }else{
        simz[s,4:6] <- c(sum$summary[1, 1], sum$summary[2, 1], sum$summary[3, 1])
      }
      
      # calculate bias
      simz$occ_bias[s] <-  (simz$occ_post[s] - occ) / simz$occ_post[s]
      simz$det_bias[s] <- (simz$det_post[s] - det) / simz$det_post[s]
      simz$hypox.p_bias[s] <- (simz$hypox.p_post[s] - hypox.p) / simz$hypox.p_post[s]
      
      # retrieve the SD 
      simz$occ_sd[s] <- sum$summary[1, 3]
      simz$det_sd[s] <- sum$summary[2, 3]
      simz$hypox.p_sd[s] <- sum$summary[3, 3]
      
      
      if (s %% 10 == 0) {
        update <- paste(Sys.time(), ' - ', dat, ' - ', s, '% done!', sep = '')
        write(update, file = 'progress.txt', append = TRUE)
      }
      
    } #end s for loop
    simz_list[[j]] <- simz 
    
  } #end j for loop
  
  save(simz_list, file= paste('~/Documents/name_of_folder/Simz_list', dat, '.Rda', sep=''))
  
} #end fxn
