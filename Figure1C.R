###Figure 1C

library(tidyverse)
library(methylKit)

load("methylKit_small.RData")

#Figure 1.2: 
c <- clusterSamples(meth.filtered, dist="euclidean", filterByQuantile =T, sd.threshold =0.75, method="ward", plot=F)

# Function to change color of the labels
labelColors = c("red","blue")
clusMember <- c("D1"=1, "D2"=1,"D3"=1,"D4"=2,"D5"=1, "D6"=2,"D8"=1,"D9"=1, "D10"=2, "D11"=1, "D12"=2,"D13"=1, "D14"=2,"D15"=2,"D16"=1, "D17"=2, "D18"=2, "D19"=1, "D20"=2) 

colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}
# using dendrapply
clusDendro = dendrapply(as.dendrogram(c), colLab)
clusDendro = dendrapply(as.dendrogram(c), colLab) %>% set("leaves_pch", 19) %>% set("leaves_col",c("red","red","blue","red","red","blue","red","blue","red","red","red","red","red","blue","blue","blue","blue","blue","blue"))

png("/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/Clustering.png", width = 300, height = 300, units='mm', res = 300)
plot(clusDendro, main = "Hierarchical Clustering", ylab = "Euclidean distance")
dev.off()