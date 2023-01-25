##Download data from https://github.com/DavideBozzi/Bozzi_et_al_2020_analysis/tree/master/R%20_analysis%20-%20Input%20Files
#Figure 1B
library(tidyverse)
install.packages("stringr")
library(stringr)

d<-read_csv("{your_dir}/Dead_fish.csv") %>% 
  dplyr::select(-Tank) %>% 
  mutate(day_num=as.numeric(str_sub(Day,-2))) 

ggplot(d, aes(x=day_num)) +
  geom_histogram()+
  geom_histogram(binwidth=1)+
  xlab("")+
  ylab("")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
