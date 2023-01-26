###Figure 2
library(methylKit)
library(tidyverse)
library(QCEWAS)
library(plotly)

load("{your_path}/MethylKit_19samples.RData")
refseq_to_ID<-read_tsv("{your_path}ID_to_refseq.tsv", col_names = c("ID","chr"))
chr_length<-read_tsv("{your_path}/GCF_000233375.1_ICSASG_v2_genomic_chrlength.txt", col_names = c("chr", "length")) %>%
  transmute(chr,cumsum=lag(c((cumsum(length))))) %>% mutate(cumsum = replace_na(cumsum, 0))

myDiff_t <- as_tibble(myDiff)
notsig<-myDiff_t %>% filter(qvalue>=0.01)
myDiff_t<-myDiff_t %>% filter(qvalue<=0.01)
myDiff.gr <- as(myDiff_t,"GRanges")
threshold_for_sig<--log10(min(notsig$pvalue))

SMP <- myDiff_t %>%
  left_join(.,chr_length, by="chr") %>%
  mutate(BPcum=start+cumsum) %>%
  
  # Calculate cumulative position of each chromosome
  dplyr::select(-cumsum) %>%
  
  # Add a cumulative position of each SNP
  mutate(hypo= ifelse((meth.diff<0),-1,1)) %>%
  mutate(SMvalue=hypo*(-log10(pvalue))) %>%
  dplyr::mutate_if(.,is.character,stringr::str_replace_all, pattern = "NW_*.*", replacement = "Scaffolds")


mydmr$DMR<-c(1:89)
DMR<-subsetByOverlaps(as(myDiff_t,"GRanges"),mydmr) %>% as_tibble()
SMP_DMR<-left_join(DMR, SMP, by=c("start"="start","meth.diff"="meth.diff")) %>%
  dplyr::mutate_if(.,is.character,stringr::str_replace_all, pattern = "NW_*.*", replacement = "Scaffolds")

m <- findOverlaps(myDiff.gr,mydmr)
gr1.matched <- myDiff.gr[queryHits(m)]
mcols(gr1.matched) <- cbind.data.frame(
mcols(gr1.matched),
mcols(mydmr[subjectHits(m)]))
gr1.matched


axisdf = SMP %>% group_by(chr) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2 ) %>% left_join(.,refseq_to_ID) %>% mutate(ID = replace_na(ID, "Scaffolds")) #mutate(ID = replace_na(ID, "Scaffolds"))

SMP$text <-paste("Chromosome: ", SMP$chr, "\nPOSITION: ",SMP$start, "\nMeth.Diff: ", SMP$meth.diff)
SMP_DMR$text <-paste("Chromosome: ", SMP_DMR$seqnames, "\nPOSITION: ",SMP_DMR$end.x, "\nMeth.Diff: ", SMP_DMR$meth.diff, "\nDMR:", gr1.matched$DMR)

p<-ggplot(SMP, aes(x=BPcum, y=SMvalue,text=text)) +
  
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.2, size=1) +
  scale_color_manual(values = c(rep(c("orange", "black"),500),"orange","orange","orange","orange","orange","orange","orange")) +
  annotate(geom = "text",label=axisdf$ID, x=axisdf$center, y=rep(c(2,0,-2),10), angle=90)+ #y=rep(c(2,0,-2),10)
  geom_point(data=SMP_DMR,aes(x=BPcum,y=SMvalue,color=as.factor(seqnames)),size=1.2,alpha=1) +
  
  
  # custom X axis:
  scale_y_continuous(labels = c(30,20,10,0,10,20,30), breaks = c(-30,-20,-10,0,10,20,30), position = "left") +     # remove space between plot area and x axis
  
  geom_hline(yintercept=threshold_for_sig, linetype="dashed", color = "red")+
  geom_hline(yintercept=-threshold_for_sig, linetype="dashed", color = "red")+
  xlab("Genomic Position") +
  ylab("-log10(p-value)") +

  # Custom the theme:
  theme_bw(base_size=16) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x  = element_blank(),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5))

ggplotly(p, tooltip="text")

P_lambda(myDiff$pvalue)
