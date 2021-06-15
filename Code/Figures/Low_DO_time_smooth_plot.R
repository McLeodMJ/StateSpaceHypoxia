DO_Data <- read.delim("~/Documents/THESIS/StateSpaceHypoxia/Data/DO_Data.txt", comment.char="#")
library(ggplot2)
library(tidyverse)

Freq_hypox <- DO_Data %>% 
  filter(DO < 1.43 & DO > 0.5) %>% 
  group_by(year) %>% 
  count()
Freq_hypox$hyp <- rep("Hypoxic", nrow(Freq_hypox))

Freq_sev_hyp <- DO_Data %>% 
  filter(DO < 0.5) %>% 
  group_by(year) %>% 
  count()
Freq_sev_hyp$hyp <- rep("Severe Hypoxia", nrow(Freq_sev_hyp))
Freq <- rbind(Freq_hypox, Freq_sev_hyp)

ggplot(Freq, aes(year, n, fill = hyp))+
  geom_smooth(model = lm)+
  labs(title="Timseries of low DO Concentrations", x="Year", y= "Count")
  #geom_col()
