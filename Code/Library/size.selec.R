# size-selextivity

size.selec <- function(dat){
  # organize minimum CW size by site
 data <- dat %>% 
   group_by(Site) %>% 
   summarise(min_CW = min(CW))
 
 # aggregate across sites
 data <- left_join(dat, data, by=c('Site'))
 data$max_CW <- rep(160, nrow(data))
 
 # eqn for linear relationship - diff btw min size obsv at 0% selec. and 1 min legal size at 100% selec.
 x <- (data$CW - data$min_CW) / (data$max_CW - data$min_CW)
 x <- ifelse(x >1, 1, x)
 #save selectivity
data$Selec.prob <- x

return(data)
 
}
