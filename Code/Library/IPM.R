
####################################
# Integral Projection Model
### evaluates size distribution across a timeseries based on SS parameters
#' fish: species - "Lingcod" "Yrock"   "Grock"   "Dsole" 
#' timestep: singular no. how long you want to run timeseries 
#' mesh= meshsize or how many bins(columns) to separate from 1 to maximum size of fish
#' N0: Stable Age distribution (SAD) 
### if 1st run --> N0= NA (or 1) - output saves nexts runs N0 as $SAD or first list [[1]]
### if 2nd run or already have inital pop. size dist. --> vector of inital sizes or [ $SAD  ]
  ### Example 1st run: N1 <- IPM() then 2nd run: N2 <- IPM( , , ,N1$SAD, )
#' cv: adds recruitment variation 
### if NO variation --> cv =0
### if including recruitment variation --> cv = number < 1

IPM <- function (fish, fi, timestep, mesh, N0, cv){

  df <- data.frame( Time = 1:timestep, 
                  Type = rep(fish, length(time)),
                  Pop.Size =  rep(NA, length(time )))


  pars <- params(fish)
  # IPM integration parameters:
    meshsize = mesh
    meshmin = 0
    meshmax = pars$Linf * 2
    
    x <- seq(from = meshmin, to = meshmax, length.out = meshsize) # different size classes
    dx = diff(x)[1] #width of the 'rectangle' to do the midpoint rule *cough cough* left rule
    
  ## size distribution of recruits
    Rvec <- dnorm(x, pars$Rlen, pars$Rlen.sd)  
    
  ## Kernel functions 
    K <- kernmat(x, pars, Dsole_f, 1)
    Fe<- fecmat(x, pars)
    
  ### Initialize the model:
    N = matrix(0, nrow = meshsize, ncol = timestep) #pop size w/ growth/mortal kernel
    E <- NULL
    Recruits <- NULL
    
   ## Prepare if SAD is known or not
    if(length(N0) == 1){
      N[,1] <- Rvec * exp(pars$R0) # SAD - not known - initialize with one pulse of new recruits
    }else{
    N[,1] <- N0 #SAD - KNOWN
    } 
    
  ### Run the model
  for (t in 2:timestep){
    N[,t] <- K %*% N[,t-1]  * dx  # midpoint rule integration
    E <- Fe %*% N[,t-1] * dx
    
  ### BH Recruitment 
    # BH eqn: (Methot and Tylor 2011 - eqn A.7)
      # Ry = 4h*R0*Eggs / S0(1-h) + Eggs(5h -1)
    Recruits <- as.vector((4 * pars$steep * exp(pars$R0) * E) / ((pars$S0 * (1 - pars$steep)) + (E * (5 * pars$steep - 1))))
    RR <- exp(rnorm(1, mean = log(Recruits) - ((cv * log(Recruits))^2)/2, sd= cv * log(Recruits) ))# change cv to 0 for no variation
    N[,t] = N[,t] + (Rvec * RR) # 
  } #end of model
  
   N0 <- N[,timestep] # save SAD
    sims <- list(SAD = N0, 
             Pop.matrix = N)
    return(sims)
}

#load("./Results/Expected_Fishing_rate.Rdata")

#source("./Library/params.R")
#source("./Library/kernmat.R")
#source("./Library/fecmat.R")
