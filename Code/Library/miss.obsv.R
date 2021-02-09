# function  to create the value 0 for lack of observation in a column 

miss.obsv <- function(dat, level, sort_col, miss_col, lev_col){

 # sp. density per trawl 
    x <- dat %>% 
      group_by(sort_col, lev_col) %>% 
      summarise(miss_col = n()) 
}
 # set species as level to report them each time
    x$lev_col <- as.ordered(x$lev_col)
    levels(x$lev_col) = level
    
 # return total combination of 4 species for each trawl ID
    all <- x %>% 
      expand(xsort_col, lev_col)
    
 # join all posibilties to df
    z <- x %>% dplyr::right_join(all)
    z$miss_col<- ifelse(is.na(z$miss_col), 0, z$miss_col) # 0 for no observation
}




#species of interest
zz <- Total_sp_data %>% 
  filter(common_name == c("Dover sole", "greenstriped rockfish", "lingcod", "yellowtail rockfish"))

# sp. density per trawl 
zzz <- zz %>% 
  group_by(trawl_id, common_name) %>% 
  summarise(density = n()) 

# set species as level to report them each time
zzz$common_name <- as.ordered(zzz$common_name)
levels(zzz$common_name) = c("Dover sole", "greenstriped rockfish", "lingcod", "yellowtail rockfish")

# return total combination of 4 species for each trawl ID
all <- zzz %>% 
  expand(trawl_id, common_name)

# join all posibilties to df
z <- zzz %>% dplyr::right_join(all)
z$Presence <- ifelse(is.na(z$density), 0, 1) #present/ absent 
z$density <- ifelse(is.na(z$density), 0, z$density) # 0 for no observation
