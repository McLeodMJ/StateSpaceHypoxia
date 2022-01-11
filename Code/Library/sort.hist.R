### START BY TRYING TO SEE WHY IT DOESNT UPDATE IN THE FOR LOOP PROPERLY 


#   Add a distribution portion
 #function (dat, trawl, mesh){
   source("./Library/params.R")
   Sp_selex <- read.csv("../Data/Sp.selex.csv") 
   load("./Results/Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors
   load("../Data/Groundfish.Hyp.Rdata")
   mesh=200
   
   # Need a vector of size bins for historgram
   pars <- params("Dsole")
   meshmax = pars$Linf * 2
   x_D <- seq(from = 1, to = meshmax, length.out = mesh)
   pars <- params("Grock")
   meshmax = pars$Linf * 2
   x_G <- seq(from = 1, to = meshmax, length.out = mesh)
   pars <- params("Yrock")
   meshmax = pars$Linf * 2
   x_Y <- seq(from = 1, to = meshmax, length.out = mesh)
   pars <- params("Lingcod")
   meshmax = pars$Linf * 2
   x_L <- seq(from = 1, to = meshmax, length.out = mesh)
   
   
   
Sp_pres <- Spe_present[Spe_present$Presence == 1, ]
Nact_pop <- matrix(NA, ncol= mesh, nrow= nrow(Sp_pres))

Nact_df <- data.frame(Trawl = Sp_pres[, 1],
                      Species = as.character(Sp_pres$common_name))
Nact_df[1:length(Nact_df), 3:202] <- NA # create space for for loop

Trawl <- NULL
g=1
y=1
l=1
d=1
for ( i in 1: nrow(Nact_df)){
if (Nact_df$trawl_id[i] == Dsole$trawl_id[d] & Nact_df$Species[i] == "Dover sole"){
   z <- Dsole[Dsole$trawl_id == Nact_df$trawl_id[i],  ]
   Trawl <- hist(z$length_cm, breaks= c(x_D, Inf), plot = F)$density # is this right because does not give a prportion by counts..
   Nact_df[i,3:ncol(Nact_df)] <- Trawl
  d = d+1
}
  if (Nact_df$trawl_id[i] == Grock$trawl_id[g] & Nact_df$Species[i] == "greenstriped rockfish"){
    z <- Grock[Grock$trawl_id == Nact_df$trawl_id[i],  ]
    Trawl <- hist(z$length_cm, breaks= c(x_G, Inf), plot = F)$density # is this right because does not give a prportion by counts..
    Nact_df[i,3:ncol(Nact_df)] <- Trawl
    g = g+1
  } 
  
  if (Nact_df$trawl_id[i] == Yrock$trawl_id[y] & Nact_df$Species[i] == "yellowtail rockfish"){
    z <- Yrock[Yrock$trawl_id == Nact_df$trawl_id[i],  ]
    Trawl <- hist(z$length_cm, breaks= c(x_D, Inf), plot = F)$density # is this right because does not give a prportion by counts..
    Nact_df[i,3:ncol(Nact_df)] <- Trawl
  y= y+1
  }
  
  if (Nact_df$trawl_id[i] == Lingcod$trawl_id[l] & Nact_df$Species[i] == "lingcod"){
    z <- Lingcod[Lingcod$trawl_id == Nact_df$trawl_id[i],  ]
    Trawl <- hist(z$length_cm, breaks= c(x_G, Inf), plot = F)$density # is this right because does not give a prportion by counts..
    Nact_df[i,3:ncol(Nact_df)] <- Trawl
  l = l+1
  }
}
