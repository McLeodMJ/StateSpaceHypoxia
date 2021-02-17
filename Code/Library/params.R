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
    R0 = 12.85,  #unfished R0 - ln(r0) or # 380777 recruits
    R0.sd = 0.35, #fixed
    steep = 0.8, 
    S0 = 469866,
    f = 0.017,
    fec.const = 0.000002805, #linear to weight
    fec.exp = 3.345,
    fec.int =1, #1
    fec.slp = 3,#0 gm-slope p# 314
    M = 0.1165,
    M.sd = 0.0056,
    K = 0.1497,
    K.sd = 0.0078,
    Linf = 47.81,
    Linf.cv = 0.114,
    Rlen = 5.405,
    Rlen.sd = 0.094 * 5.405, 
    Mat.len = 35,  #p. 34
    Mat.slp = -0.775)# page 314
    
  }
  
  # Parameters from Northern model & Female
 if(sp == "Lingcod"){
    Params <- list( 
      R0 = 9.0669, #ln(r0) // #8037
      R0.sd = 0.16454, #have 95% CI
      steep = 0.7,
      S0 = 37974, #unit is mt
      f = 0.115,
      fec.const = 0.00000276,
      fec.exp = 3.28,
      fec.int =1,
      fec.slp = 3,
      M = 0.257,
      M.sd = 0.4384,
      K = 0.0173,
      K.sd = NA,
      Linf = 108.6,
      Linf.cv = 0.06,
      Rlen = 17.2792,
      Rlen.sd = 0.1436 * 17.28, #cv-length min * mean #page 97
      Mat.len = 56.7, # p. 97 but 64 on #page 30
      Mat.slp = -0.269) #p.97
    
  }
  
  # Female parameters
 if(sp == "Grock"){
    Params <- list( 
      R0 = 9.62,# - ln(R0) // 15041 recruits
      R0.sd = 0.15, # 95 CI: 2073-109,131
      steep = 0.69,
      S0 = 7090,
      f = 0.0682, #exploitation rate: SPR(msy)
      fec.const = 371200, #fec. intercept
      fec.exp = 63300, #fec. slope
      fec.int =1, #1 
      fec.slp = 3, #0 p#207
      M = 0.08, 
      M.sd = NA, # 95 AI: 0.02-.4
      K = 0.11,
      K.sd = 0.003,
      Linf = 33.67,
      Linf.cv = 0.07,
      Rlen = 9.25,#page 74
      Rlen.sd = 0.09 * 9.25, #cv given
      Mat.len = 20.97 ,  #p. 207
      Mat.slp = -0.66) #p. 207
    
  }
  
  # Parameters from Northern model & Female
 if(sp == "Yrock"){
    Params <- list( 
      R0 = 10.83, # ln(R0) // 49090000 recruits
      R0.sd = 0.15, #95 ci: 17.86-134.94 mil
      steep = 0.718,
      S0 = 14.996, #trillion eggs
      f = 0.089, #95 CI: 0.085-0.093
      fec.const = 1.1185e-11,
      fec.exp = 4.59,
      fec.int =1, #0
      fec.slp = 3, #3 p# 80
      M = 0.173,
      M.sd = NA,
      K = 0.17,
      K.sd = NA, #literature values used/ do not give SD
      Linf = 52.2,
      Linf.cv = 0.1, #do not give SD or CV - SD estimated at 0.1
      Rlen = 14.689,#page 80
      Rlen.sd = 0.105 * 14.689,# given as cv
      Mat.len = 42.49,  #p. 80
      Mat.slp = -0.401) #guess
   }
  return(Params)
}
