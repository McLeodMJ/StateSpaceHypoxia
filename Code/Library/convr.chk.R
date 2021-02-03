########################################################################################
### Function to exclude posteriors with convergence 
# summ = sum$summary - the summary table from MCMC that is the average across chains
# x = empty vector of zeros to use as an indicator of convergence 
## - if x has any NA's then convergence and do not report values 
#########################################################################################

convr.chk <- function(summ, x){
  
   for(n in 1: nrow(summ)){
      # insufficient rhat >=1.05
      if(summ[n, 10] >= 1.05){
        x[n] = NA
        x <- x[is.na(x)]
      } 
     
      #insufficient effect size <= 100
      if(summ[n, 9] <= 100){
        x[n] = NA
        x <- x[is.na(x)]
      }
    } 
    return(x)
}
