Stan!
  

library(rstan)
rstan_options(auto_write = TRUE) #tany suggestions made from library rstan will be added auto. to help smooth the process

- store stan model as text file to call model 

- data is stored in list
list(N= 100, X=X) #stores number of observations 1st then columns in df as list 

#Want to save the stan.model
save(model, file = "file.name")
  - will NOT be able to look at statistics after the model has been saved 

load('file.name') #reloads the model to environment  

each line of code in stan ends with a semicolon ;
------------------------------
  
Prob. 

- conditional prob. - a | b OR a given b
- bayes theorm 
posterior: what you are interested in finding (b|a)
  -p(b|a)
  - b= nodel param. a= data
likelihood: trying the possible values for model param 
  -p(a|b)
prior: infor that you have about the param before collecting data
  -p(b)
marginal posterior (or evidence): integrates all possible values of param. b.
  -p(a)
  - makes sure that everything sums to one so you get a prob.#
  - normalization term SO we can pften exclude EXCEPT for model comparison 


---------------

Stan model
###### three main blocks

data {
  1. declare the no. of observations (n)
  2. declare the data (arrays of size n)
    - include any uncertainty = sigma
}
parameters{
  - declare parameters (unobserved variables)
  - prob. need to have prob boundaries 0-1
  - mu = population tratment effect #when hierarchal models
    - must scale the population 
}
model{ 
  1. theta ~ MEANS theta is distributed from .... or theta is drawn from.....
  - where you declare the distrubtions pirors
  2. define the likelihoods
    - if present absent or clicked or didtnt click (Ctr) then bernoulli(theta)
    - y1 ~ dist() of likelihood
}

# USe
print() # to make sure everything is syntactically correct --> compile 
# EX like fxn where you should alwats do it at the end to check 

#OPTIONAL BUT wil help later on?
generated quantities {
  #variables interested in from previously declared variables 
  real Delta_theta; 
  Delta_theta = theta1 - theta2 
}

--------------------
  
 NEXT fit data to the model 
1. turn data into list of each column being a list
  -EX. data <- list(x,y, n = length(n.obs))
2. fit <- stan(file = 'stan.model', data=data)
  - extract parameters from 'fit'
    - params = extract(fit)
3. observe dist. 
  - plot(density(params$theta1))
    - if gap of of peak is large then large uncertainty 
  - plot(density(params$Delat-theta)) #look at the difference of comparable params
    - 0 on x means that the params are the same but...
    if it shifts to left or right then that param it is shifted to is better
-------------------------------

CONVERGENCE!
  - CLT: if sample is sufficiently large, random and independent then the dist of sample means is approx. normal
  - MCMC method used to draw samples frrom post dist and useful when nontrivial to pull indep. 
    - assumes CLT= MCMC is a normal dist --> convergence
    - convergence: reached target dist. and not biased or correlated 
  
- run the stan.model for multiple chains
    - each chain is a different dist. and has residual effects
    - discard initital iterations = 'burn-in'
  - sample pass of chains can be3 visually inspected 
    -  traceplot(fit, pars= c('mu', 'tau'))
      - change the fit variable by running the model with more iteration sto improve mixing 
        - too many parameters would be challenging to do traceplots for each so.....

  r-hat
    - R = v/w 
      - v= pooled variance of all chains 
      - w = within chain variance 
  - - if samples mixed well will be == 1
      - WANT r.hat 1- 1.1
      - rhat > 1.1 = varaiance of the combined change will be grEater than thevariance ind. chains
- print (fit) #gives the different rhats

--------------------------------
    
effective sample size 
  - the no. of indep. samples as determined from the increase in uncertainty in posterior est. from correlated samples 
    - summarh(fit)$summary[, 'n_eff'] --> ratio of effective samples / total iterations 
      - want the ratio to be > 0.01 otherwise has bias
  

ERROR: bulk_ess or tail_ess
- tails of the posterior
  - monitor(extract(fit, permuted= false, inc_warmup=FALSE))
    - having longer iterations should help make the errors disapear


Divergence
- CALC. THE DERIVATIVES OF THE POSTERIR DIST FXN 
  - HAMILTONIAN DYNAMICS 
  - divergences errors are trying to tell you that it cannot plot along the curve of tbe dist
    - pairs(fit, pars=c("mu", "tau", 'lp__'))
  - use a smaller step size to geta more accurate trajectory 
    - fit = stan('model.stan', data= data, chains =4, iter= 1000, warmup= 500, control= list(adapt_delta=0.85))
      - if the divergence does NOT go down with ^ step sizes then there is something wrong with priors or need to reparameterize


-----------------------------
  
Homogenity: subsamples of datasets have the same statistical properties as full datasets
  - homosecadacity data = uniformly dist. (same variance)
    - evenly distributed around x
heterogenity: different stat. prop. of full dataset
  - variance of Y ^^ as a fuction of covariate X
  - combining datasets --> relationship may vary subsamples (discrete hetero.)
    - use hierarchical model

- homoscedastic model 
  - data = n.obs & variates and covariates 
  - param = params trying to define (i.e. alpha, beta && sigma)'
'
  - model = priors (dists)
  - likelihood = variate ~ covariates, sigma #ASSUME constant variance 

- heteroscedastic model 
- data = n.obs & variates and covariates 
- param = params trying to define (i.e. alpha, beta && sigma)
- model = priors (dists)
- likelihood = variate ~ covariates, sigma*X #NOT constant variance
  - adds dependence of covariance 

## use pairs to check 

- time series data = non-normal heteroscedastictiy where scale of noise varies over time 

- modeling heterosc.
  1. plot percent change over time
    - IF heterogen (step 2)
  2. 

  ----------------------------------------------
  ##  FABIO NOTES 
    priors for the optim section 
  
  - inits: gets all the priors and draws on them through dist --> feeds the model 
  
  - vbgm( call parameters for all linf, tnaught, k) --> then make the vbgm into a list where separate 
  
  CONSTRAINTS 
  - real for estimated singular values 
  - int is for singular numbers based on data 
  - vector is a list of numbers 
  - <lower=0> constrains value to positive number
  
  need real in the parameters 
  - need to include param. est in the data section 


  -------------------------------------------
    # CHAPTER NOTES
require(rethinking)
   - ulam() fxn --> translates the data list to get samples from the post. --> and converts to stan model for you 
  
WARM UPS / ITERATIONS
- Not the same as burn ins b/c will heat things up then once real iterations then should go straight to target dist. 
 1. must figure out optimal effective number of samples 
 2. look at the est. of of n_eff from the post. dist. 
      - 
 3. if you want more than posterior means then need more samples to get exact shape of post. dist.
 4. STAN error: tail ESS - means that the quality of intervals to be questionable 
        - DO more samples
 *** generally rule: if havingn trouble w/ getting post. dist. ^^^ samples, if not trouble, reduce it the # of samples
  
 
CHAINS
 1. start with one chain - some errors only show up with one chain
 2. traceplot & trankplot to see that the chain works 
 3. if chains are valid --> need replication with more than 1 chain
 4. make inferences from chain or chains
  
  
** 1 short chain to debug --> 4 chains for verification and inference 

rhat: measure of convergence 
  - want rhat slightly greater than 1 [1-1.1]
    - if CRAZY rhat then run a pairs plot to see issue or traceplot
    - change to weaker regularizing priors if NOT mixing or stationary 
      - EX. instead of dnorm(0,1000) do dnorm(0,10)

  
Error: treedepth - indicative of inefficient chains 
  - control=list(max_treedepth=15)





####################### results #######################################
# look at diagnostics 
rstanarm::launch_Shinystan(model)

#get summarys via R
summary(model)$summary %>% head()

#plots of param.
plot(model)

#Divergence check
rstan::check_divergences(model)

#looking at results per chain
rstan::get_sampler_params(model) %>% #returns a list of one object per chain
  set_names(1:n_chains) %>% #sets it so that name of chains is 1-> number of them
  map_df(as_data_frame,.id = 'chain') %>% 
  group_by(chain) %>% 
  mutate(iteration = 1:length(chain)) %>% 
  mutate(warmup = iteration <= warmups) # OR inc_warmup=FALSE to remove warm ups

# want the HMC to test new parameters on reasonable timeline 
max_treedepth = 10 #Default
# run a ggplot to see if the parameters are being checked only ecause it hit maxtree_depth each time --> New higher value 

#stepsize = resolution 
ggplot(model, aes(iteration, stepsize__, color = chain)) + 
  geom_line()  + 
  scale_color_locuszoom()

############ still diverge4nce?

1. feature engineering 
  - when each param that we are looking at for dependence are on completely different numerial scalesss
    - EX> 




  
  