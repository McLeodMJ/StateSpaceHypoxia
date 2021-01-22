#################################################################
### Function for storing and calling species specific parameters
# sp = name of species: "Dsole", 
#'                      "Lingcod",
#'                       "Yrock", 
#'                       "Grock"
#################################################################

params <- function(sp){
  # Female Parameters
  if(sp == "Dsole"){
   Params <- list(
    rec = 12.85,  #unfished R0 - ln(r0) or # 380777 recruits
    rec.sd = 0.35, #fixed
    steep = 0.8, 
    s0 = 469,866,
    sy = 393,507,
    f = 0.017,
    fec.const = 0.000002805, #linear to weight
    fec.exp = 3.345,
    m = 0.1165,
    m.sd = 0.0056,
    k = 0.1497,
    k.sd = 0.0078,
    Linf = 47.81,
    Linf.sd = (0.114 * 47.81) )#sd = .114
    
  }
  
  # Parameters from Northern model & Female
 if(sp == "Lingcod"){
    Params <- list( 
      rec = 9.0669, #ln(r0) // #8037
      rec.sd = 0.16454, #have 95% CI
      steep = 0.7,
      s0 = 37974, #unit is mt
      sy = 21974, # unit not  reported by probs mt
      f = 0.115,
      fec.const = 0.00000276,
      fec.exp = 3.28,
      m = 0.257,
      m.sd = 0.4384,
      k = 0.0173,
      k.sd = NA,
      Linf = 108.6,
      Linf.sd = (0.06 * 108.6)) #cv= .06
    
  }
  
  # Female parameters
 if(sp == "Grock"){
    Params <- list( 
      rec = 9.62,# - ln(R0) // 15041 recruits
      rec.sd = 0.15, # 95 CI: 2073-109,131
      steep = 0.69,
      s0 = 7090,
      sy =2836,
      f = 0.0682, #exploitation rate: SPR(msy)
      fec.const = 371200, #fec. intercept
      fec.exp = 63300, #fec. slope
      m = 0.08, 
      m.sd = NA, # 95 AI: 0.02-.4
      k = 0.11,
      k.sd = 0.003,
      Linf = 33.67,
      Linf.sd = (0.07 * 33.67)) #cv is 0.07
    
  }
  
  # Parameters from Northern model & Female
 if(sp == "Yrock"){
    Params <- list( 
      rec = 10.83, # ln(R0) // 49090000 recruits
      rec.sd = 0.15, #95 ci: 17.86-134.94 mil
      steep = 0.718,
      s0 = 14.996, #trillion eggs
      sy = 11.278, #trillion eggs
      f = 0.089, #95 CI: 0.085-0.093
      fec.const = 1.1185e-11,
      fec.exp = 4.59,
      m = 0.173,
      m.sd = NA,
      k = 0.17,
      k.sd = NA, #literature values used/ do not give SD
      Linf = 52.2,
      Linf.sd = 3) #do not give SD - estimated from figure 5 (northern/ female)
   }
  return(Params)
}
