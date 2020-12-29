setwd("~/Box Sync/McLeod_thesis/code")
DOdata <- read.delim("~/Box Sync/McLeod_thesis/Data/DO_Data.txt", comment.char="#")
library(tidyverse)
library(pracma)

plot(DO~time, data= DOdata[DOdata$year == 2020,], col= "red", type='l')
abline(h=1.43)

data <- as_tibble(DOdata)
Low_hypox <- data %>%
  select(DO, year) %>%
  filter(DO <= 1.43) %>%
  group_by(year) %>%
  tally() #low years: 2018 and 2020

High_Hypox <- data %>%
  select(DO, year) %>%
  group_by(year) %>%
  mutate(mean = mean(na.omit(DO))) %>% 
  tally()
 #high years 2010 & 2014


ggplot(DOdata, aes(x= time, y=DO, color=year))+
  geom_line()+
  scale_fill_brewer(palette = "Accent")+
  facet_wrap(~year, 4)
#2016 * 2020 start around 150-200 time (not complete series)


# make X into a timeseries object
Dt <- ts(data=DOdata$DO, start = min(DOdata$year), end = max(DOdata$year)) # note that this assumes sequential evenly-spaced sampling. If that is not the case we have to do some padding...

# scale to unit variance
Dt <- Dt/sd(Dt, na.rm=TRUE)

# detrend
Time = min(DOdata$year) : max(DOdata$year)
m <- lm(Dt~Time)
Dt2 = m$resid
Dt2 = Dt2/sd(Dt2) # rescale to unit variance

# Now do the FFT
f <- seq(from=0, to=0.5, length.out = ceiling(length(Dt2) / 2)) # vector of frequencies
F <- abs((fft(Dt2) ) ) / length(Dt2)
F = F[1:length(f)]

# make into dataframe for ggplot
DsD <-  data.frame(x=1/f, y= F^2) 
#plot(1/f, F, type='l', xlim=c(0,20))


# get ar scale
Dt.ar <- ar(Dt2, na.action=na.pass)
Dt.ar #0.4952


# Red noise significance cutoff (based on Torrence & Compo 1998)
lag1 = 0.4952  # Get this value from the ar() operation above #coefficient correlation btw successive pts
fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(f*2*pi) + lag1^2) #eqn 16 w/in cos() should have / length(f) # No - f is already equal to k/N in their notation
fft_theor <- var(Dt2, na.rm=TRUE) * fft_theor / (2 * length(f)) #rescaling for unit variance?
chi2 = qchisq(0.05, 2, lower.tail=FALSE) / 2
signif = chi2*fft_theor

DsD$signif = signif

ggplot(DsD, aes(x,y))+
  geom_line()+
  geom_line(aes(x,signif, color= "red"))

# Simulating new datasets based on the FFT

# We do this by randomizing the phase (imaginary part) of the FFT spectrum
Theta <- runif(n=length(f)) * 2 * pi # 2pi converts to radians
# Now do inverse FFT with modulus (real part) of original FFT but randomized phase
Z = Re(fft(Dt2[1:length(f)])) * exp(1i*Theta) #DIFFERENT LENTHS BUT TOOK FIRST 6 LIKE ON LINE 28
ZZ = Re(ifft(Z))
# repeat as many times as you like...


# An example just to demonstrate...
t = 1:1000
x1 = sin(2*pi*t/10)
x2 = sin(2*pi*t/100+5)
x3 = runif(length(t))
x = x1 + 0.5*x2 + x3 #equivalent to m ?
plot(t,x,type='l')

lag1 = 0.23  # Get this value from the ar() operation above #coefficient correlation btw successive pts
fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(f*2*pi) + lag1^2) #eqn 16 w/in cos() should have / length(f) # No - f is already equal to k/N in their notation
fft_theor <- var(Dt2, na.rm=TRUE) * fft_theor / (2 * length(f)) #rescaling for unit variance?
chi2 = qchisq(0.05, 2, lower.tail=FALSE) / 2
signif = chi2*fft_theor


Theta = runif(n=length(t))*2*pi
Z = Re(fft(x))*exp(1i*Theta)
ZZ = Re(ifft(Z))

plot(ZZ,type='l')

# check the spectrum
f = seq(from=0,to=0.5,length.out=ceiling(length(ZZ)/2)) # vector of frequencies
F0 = abs(Re(fft(x)))/length(x)
F0 = F0[1:length(f)]
F = abs((fft(ZZ)))/length(ZZ)
F = F[1:length(f)]

plot(1/f,2*F^2,type='l',xlim=c(0,100))
lines(1/f,2*F0^2,type='l',col='red')




