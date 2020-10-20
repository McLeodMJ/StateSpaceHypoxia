# Hypoxia STAN Model

##########################################
#psi = prob. space is occupied by thst sp.
#p1 = prob. of detection w/ covariate
#hypox.p = variation around hypoxia
##########################################
setwd("~/Documents/THESIS/StateSpaceHypoxia/Code")
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
load("../Data/Random_hypoxia.Rdata") #remove /.Data/


# Generate a random simulation 
n= 1000 #no. of simulations  
sim <- data.frame( trawl = seq(1, n, by=1), 
                   pres = rep(NA, n), 
                   hypox = rand) #pulls a random sample from WCGBTS - Needs to be improved to incorperate a fourier trans.

#Severe hypoxia as logistic regression w/ 0.965 as midpoint of hypoxia 
sim$p <- 1/ (1+exp(-(sim$hypox - 0.965)))
plot(p~hypox, sim )

sim$pres <- rbinom(nrow(sim), 1, sim$p) #random presence/absence

# Functions for transformations
logit <- function(x){
  log(x/(1-x))
}
inv_logit <- function(x){
  exp(x)/(exp(x)+1)
}

# Initial variables 
N = 1000
psi = 0.5
hypox = sim$hypox 
p = 0.2
hypox.p = 1 # not sure what to do for this value which is used below in the inv-logit eqn 

# factoring in hypoxia as a covariate 
p.full <- inv_logit(logit(p) + hypox.p * hypox)
# encounters [0 or 1]
enc <- rbinom(N,1,0.5) * rbinom(N, 1, p.full)


### STAN MODEL
occ.mod <- "
data{
	int<lower=1> N; //number of samples
	int<lower=0> enc[N]; //encounters
	vector[N] hypox; //hypoxia
}
parameters{
	real logit_psi; //logit occupancy
	real logit_p; //logit detection
	real hypox_p; //hypoxia slope
}
transformed parameters{
	real<lower=0, upper=1> psi;
	real<lower=0, upper=1> p;

	psi = inv_logit(logit_psi);
	p = inv_logit(logit_p);
}
model{
  // uninformative priors
	logit_psi ~ normal(0,10);
	logit_p ~ normal(0,10);
	hypox_p ~ normal(0,10);


	for(i in 1:N){
		if(enc[i] > 0){
			//the site was occupied, and you detected it
			target += log_inv_logit(logit_psi) + bernoulli_logit_lpmf(1| logit_p + hypox_p*hypox[i]);
		}else{
			//the site was occupied but you didn't detect it & the site was unoccupied
			target += log_sum_exp(log_inv_logit(logit_psi) + bernoulli_logit_lpmf(0| logit_p + hypox_p*hypox[i]), log1m_inv_logit(logit_psi));
		}
	}
}"

# Creating data in list format necessary for Stan 
occ.stan.dat <- list(N = N,
                     enc = enc,
                     hypox = hypox)

# Run Stan model 
occ.stan <- stan(model_code = occ.mod,
                 pars = c('psi','p','hypox_p'),
                 data = occ.stan.dat,
                 chains = 1,
                 warmup = 500,
                 iter = 1000,
                 #init = inits, #may need to be species specific 
                 control = list(adapt_delta = 0.95))
                 
# Check Stan model runs
pairs(occ.stan)
traceplot(occ.stan)

