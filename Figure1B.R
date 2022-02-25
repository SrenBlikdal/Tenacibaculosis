###Figure 1B

library(tidyverse)
library(methylKit)
library(ggrepel)

load("methylKit_small.RData")

sample.id <- c("D1","D2","D3","D4","D5","D6","D8","D9","D10", "D11","D12","D13","D14","D15","D16","D17","D18","D19","D20")
treatment_n <- c(1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0)
treatment <- gsub(1,"Sick",treatment_n) %>% gsub("0","Healthy",.)
design.matrix <- bind_cols(ID=sample.id, treatment=treatment)

pca <- percMethylation(meth.filtered) %>% t(.) %>% prcomp(.,center=T, scale. = F)
percent_variance_pca<-summary(pca)$importance["Proportion of Variance",] * 100

pca_plot<-as_tibble(pca$x) %>% bind_cols(design.matrix)

png(file="/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/PCA.png", width = 300, height = 600, units='mm', res = 300)
ggplot(pca_plot,aes(x=PC1, y=PC2, color=treatment,label= ID)) +
  geom_point(size=3)  +
  geom_text_repel(size=7) +
  scale_color_manual(values=c("blue", "red"))+
  xlab(label = paste("PC1 (",round(percent_variance_pca[1],1),"percent of the total variance)")) +
  ylab(label = paste("PC2 (",round(percent_variance_pca[2],1),"percent of the total variance)"))+
  ggtitle("PCA \n 5X threshold")+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
  theme_bw(base_size=16)+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend("Disease status"))
dev.off()