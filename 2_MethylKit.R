#!/usr/bin/Rscript

#load the needed libraries
library(tidyverse)
#tidyverse_1.3.0 
library(sessioninfo)
#sessioninfo_1.1.1
library(methylKit)
#methylKit_1.16.1 
library(edmr)
#edmr_0.6.4.1

setwd("{Dir_with_bismark_cov_files}")

file.list<-list(
"D1.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D3.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D5.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D6.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D8.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D9.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D10.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D11.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D12.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D13.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D14.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D15.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D16.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D17.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D18.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D19.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D20.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz") 

my_obj<-methRead(file.list, sample.id = list("D1","D2","D3","D4","D5","D6","D8","D9","D10", "D11","D12","D13","D14","D15","D16","D17","D18","D19","D20"), assembly = "ICSASG_v2", pipeline = "bismarkCoverage", treatment= c(1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0), context= "CpG", mincov=5)
#treatment is meta data from Bozzi et al records:0 is healthy/Mycoplasma, 1 correspond to sick/Avivibro. 
#context =CpG, means we only look at CpG sites 

filtered.my_obj<-filterByCoverage(my_obj,lo.count=5,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)

meth.filtered<-unite(filtered.my_obj, destrand=F)

myDiff<-calculateDiffMeth(meth.filtered, mc.cores=40)

mydmr<-edmr(myDiff, dist=150,DMR.methdiff = 0,DMC.methdiff = 0, num.DMCs = 10)

#Remove intermediate files 

rm(my_obj)
rm(filtered.my_obj)

#Save workspace
save.image(file ="{outputdir/MethylKit_19samples.RData}")

q()
n
