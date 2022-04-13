require(ggplot2)

DOdata <- read.delim("~/Box Sync/McLeod_thesis/Data/DO_Data.txt", comment.char="#")
load("DO_sims.Rdata")



ggplot(dat_NA[dat_NA$year %in% 2010:2015, ], aes(Jday, DO, color=year))+
  geom_line()+
  scale_fill_brewer(palette = "Accent")+
  facet_wrap(~year, 4)+
  geom_hline(yintercept = 1.43, color= "red")+
  ggtitle("Summertime DO at a bottom mooring off Newport")+
  xlab("Julian Days")+
  ylab("DO Concentrations")+
  guides(color="none")


ggplot(DsD, aes(Period, Freq_sq))+
  geom_line()+
  geom_line(aes(Period, Signf, color= "red"))+
  facet_wrap(~Year)+
  ggtitle('Fourier spectra of DO timeseries, used to generate realistic simulated timeseries for method development')+
  guides(color="none")
