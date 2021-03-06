###Figure 1A

library(tidyverse)
library(methylKit)
library(GenomicRanges)

load("methylKit_small.RData")

mean.meth_H<-reorganize(meth.filtered,sample.ids=c("D4","D6","D10","D12","D14","D15","D17","D18","D20"), treatment=c(0,0,0,0,0,0,0,0,0)) %>% percMethylation(.) %>% rowMeans(.)
mean.meth_S<-reorganize(meth.filtered,sample.ids=c("D1","D2","D3","D5","D8","D9","D11","D13","D16","D19"), treatment=c(1,1,1,1,1,1,1,1,1,1)) %>% percMethylation(.) %>% rowMeans(.)

mean.meth_H_GR<-bind_cols("chrom"=meth.filtered$chr,"chromStart"=meth.filtered$start,"chromEnd"=meth.filtered$end, "mean.meth"=mean.meth_H) %>% as(.,"GRanges")
mean.meth_S_GR<-bind_cols("chrom"=meth.filtered$chr,"chromStart"=meth.filtered$start,"chromEnd"=meth.filtered$end, "mean.meth"=mean.meth_S) %>% as(.,"GRanges")

inputnames<-c(seq(-6000, -50, by = 50),seq(50, 4000, by= 50))

headder<-c("chrom","chromStart","chromEnd","name", "score", "strand")

H<-inputnames %>% map(function(x) read_tsv(paste0("~/Documents/Master_thesis/annotation/GCF_000233375.1_ICSASG_v2_genomic_79000_",x,"promoter.bed"), col_names = headder) %>% as(.,"GRanges") %>% subsetByOverlaps(mean.meth_H_GR,.) %>% as_tibble() %>% transmute(mean=mean.meth, pos=x))
S<-inputnames %>% map(function(x) read_tsv(paste0("~/Documents/Master_thesis/annotation/GCF_000233375.1_ICSASG_v2_genomic_79000_",x,"promoter.bed"), col_names = headder) %>% as(.,"GRanges") %>% subsetByOverlaps(mean.meth_S_GR,.) %>% as_tibble() %>% transmute(mean=mean.meth, pos=x))

all_H<-bind_rows(H) %>% group_by(pos) %>% summarise(., "m"=mean(mean), "sd" = sd(mean)) %>% mutate(., DS = "Healthy")
all_S<-bind_rows(S) %>% group_by(pos) %>% summarise(., "m"=mean(mean), "sd" = sd(mean)) %>% mutate(., DS = "Sick")

p<-bind_rows(all_H,all_S)

#png("/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/mean_meth.png", width = 1000, height = 400, units='mm', res = 300)
ggplot(p, aes(x=as.numeric(pos),y=m, color=DS))+
  geom_line(size=2)+
  ylab("Mean methylation")+
  xlab("Position from TSS")+ 
  scale_color_manual(values=c("blue", "red"))+
  theme_bw(base_size=16)+
  guides(color=guide_legend("Disease status"))
#dev.off()