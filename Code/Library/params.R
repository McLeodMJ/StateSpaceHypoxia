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
     name = "Dsole",
    R0 = 12.85,  #unfished R0 - ln(r0) or # 380777 recruits
    R0.sd = 0.35, #fixed
    steep = 0.8, 
    S0 = 469866,
    f = Sp_depl[Sp_depl$Species == "Dsole",3], #hybrid page 315 -0.04
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
    Mat.slp = -0.775, # page 314
    selc.50 = 35,  # pg.35 length at 50% mature
    selc.100 = 50, # page 100 figure 26
    depl = 0.837,#page 7
    selc.p1 = 33.120,  # Dbl Norm. table 21 p. 66 [peak at .6 then odd decline]
    selc.p2 = 1.697, 
    selc.p3 = 3.928, 
    selc.p4 = 0.665,
    selc.p5 = NA,
    selc.p6 = 1.844,
    selc.F.p4 = -0.723,
    selc.F.p6 = -1.455,
    bin.len = 65,# page 37 & 207 of assessment
    bin.wth = 1,
    selc = spline(Sp_selex$Dsole, n=200)$y)
    
  }
  
  # Parameters from Northern model & Female
 if(sp == "Lingcod"){
    Params <- list( 
      name = "Lingcod",
      R0 = 9.0669, #ln(r0) // #8037
      R0.sd = 0.16454, #have 95% CI
      steep = 0.7,
      S0 = 37974, #unit is mt
      f = Sp_depl[Sp_depl$Species == "Lingcod",3] , #webplot digit page 238 -0.117
      fec.const = 0.00000276,
      fec.exp = 3.28,
      fec.int =1,
      fec.slp = 3,
      M = 0.257,
      M.sd = 0.4384,
      K = 0.1283,#page 97
      K.sd = NA,
      Linf = 108.6,
      Linf.cv = 0.06,
      Rlen = 17.2792,
      Rlen.sd = 0.1436 * 17.28, #cv-length min * mean #page 97
      Mat.len = 56.7, # p. 97 but 64 on #page 30
      Mat.slp = -0.269, #p.97
      selc.50 = 56.63,  # p. 40
      selc.100 = 140, # pg 233 figure 110
      depl = 0.579, #page 8 in decimal of %
      selc.p1 = 61.2144,  # Dbl Norm. table 7 p. 99 [dble-norm]
      selc.p2 = -15, 
      selc.p3 = 6.45783, 
      selc.p4 = 7.05119, 
      selc.p5 = -999, 
      selc.p6 = -999,
      bin.len =140, #data.ss_new file from SS
      bin.wth = 2,
      selc = spline(Sp_selex$Lingcod, n=200)$y)
  }
  
  # Female parameters
 if(sp == "Grock"){
    Params <- list( 
      name = "Grock",
      R0 = 9.62,# - ln(R0) // 15041 recruits
      R0.sd = 0.15, # 95 CI: 2073-109,131
      steep = 0.69,
      S0 = 37947,
      f = Sp_depl[Sp_depl$Species == "Grock",3] , #pope page 228 -0.1
      fec.const = 371200, #fec. intercept PAGE 17 & 54
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
      Mat.slp = -0.66, #p. 207
      selc.50 = 21,# p. 68 L50%
      selc.100 = 28, # p. 68 L100%
      depl = 0.814, #page 8 in decimal of %
      selc.p1 = 26.77,  # Dbl Norm. table 33 p. 74
      selc.p2 = NA, 
      selc.p3 = 4.03, 
      selc.p4 = NA,
      selc.p5 = -4.68,
      selc.p6 = NA,
      
      selc.p1 = 21,  # p 214 the table of params [init] - asymp.
      selc.p2 = 3, 
      selc.p3 = 3.7, 
      selc.p4 = 6,
      selc.p5 = -5,
      selc.p6 = -999,
      bin.len = 45, #max length page 145
      bin.wth = 1,
      selc = spline(Sp_selex$Grock, n=200)$y) 
    
  }
  
  # Parameters from Northern model & Female
 if(sp == "Yrock"){
    Params <- list( 
      name = "Yrock",
      R0 = 10.83, # ln(R0) // 49090000 recruits
      R0.sd = 0.15, #95 ci: 17.86-134.94 mil
      steep = 0.718,
      S0 = 14.996, #trillion eggs
      f = Sp_depl[Sp_depl$Species == "Yrock",3] , # page 22 expl. rate - 0.01169 [ Andi Stephens confirmed 2/25/21 that it is pope's eqn & recorde as exp.]
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
      Rlen = 14.689, #page 80
      Rlen.sd = 0.105 * 14.689, # given as cv
      Mat.len = 42.49,  #p. 80
      Mat.slp = -0.401,
      selc.50 = 42.49, # page 30
      selc.100 = 53.58, #page 80
      depl = 0.752, #page 8 in decimal of %
      selc.p1 = 50.169,  # Dbl Norm.eqn table 9 p. 82 (asympt)
      selc.p2 = 70, 
      selc.p3 = 4.541, 
      selc.p4 = 70, 
      selc.p5 = -999, 
      selc.p6 = -999 ,
      bin.len = 65, # found in YTRK.north.data.ss
      bin.wth = 1,
      selc = spline(Sp_selex$Yrock, n=200)$y)
   }
  return(Params)
}



# north model for lingcod and yellowtail 
# all female parameters
