###Figure0
library(tidyverse)
OTU<-read_csv("~/Documents/Master_thesis/Crappyfish/Crappyfish_OTU.csv")

OTU_MB<- OTU %>% filter(Tissue== "M", Time== "B") %>% 
  rowid_to_column(., "EpiID") %>%
  arrange((`OTU2(Mycoplasmataceae)`)) %>% 
  rowid_to_column(., "ID") %>%
  mutate(other=10000-(`OTU1(Alivibrio)`+`OTU2(Mycoplasmataceae)`)) %>% 
  dplyr::select(.,EpiID,ID,Health, `OTU1(Alivibrio)`, `OTU2(Mycoplasmataceae)`,other) %>%
  gather(key="OTU", value="number", `OTU1(Alivibrio)`, `OTU2(Mycoplasmataceae)`,other)

OTU_MB %>%
  filter(Health=="H")%>%
  ggplot(.,aes(x= as.factor(ID),y=number,label=Health, fill=factor(OTU,levels = c("other","OTU1(Alivibrio)","OTU2(Mycoplasmataceae)")), width=0.7)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual("OTU", values = c("OTU1(Alivibrio)" = "orange2", "OTU2(Mycoplasmataceae)" = "skyblue2", "other" = "gray43"))+
  xlab("Disease status")+
  ylab("Relative aboundance")+
  coord_flip()+
  ggtitle("16S metabarcoding  \n N=40")+
  theme(panel.background = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y.left =element_blank(),
        axis.ticks.y=element_blank(), 
        axis.text.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

