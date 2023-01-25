###Figure 2C

library(tidyverse)
library(methylKit)

load("{your_path}/MethylKit_19samples.RData")

c <- clusterSamples(meth.filtered, dist="euclidean", filterByQuantile =T, sd.threshold =0.75, method="ward", plot=F)

# Function to change color of the labels
labelColors = c("red","blue")
clusMember <- c("D1"=1, "D2"=1,"D3"=1,"D4"=2,"D5"=1, "D6"=2,"D8"=1,"D9"=1, "D10"=2, "D11"=1, "D12"=2,"D13"=1, "D14"=2,"D15"=2,"D16"=1, "D17"=2, "D18"=2, "D19"=1, "D20"=2) 

colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = labCol, col=labCol, pch = 19,cex=1.5))
  }
  n
}
# using dendrapply
clusDendro = dendrapply(as.dendrogram(c), colLab)

plot(clusDendro, main = "Hierarchical Clustering", ylab = "Euclidean distance")

