library(tidyverse)

load("~/Documents/Crappyfish/MethylKit_19samples.RData")

### Run in R
export.bed(mydmr, con = "~/Documents/Crappyfish/dmrs.bed")

##This part I run in the terminal with input files from Reference_annotation 

for feature in GCF_000233375.1_ICSASG_v2_genomic_250promoter.bed GCF_000233375.1_ICSASG_v2_genomic_750promoter.bed GCF_000233375.1_ICSASG_v2_genomic_5000promoter.bed GCF_000233375.1_ICSASG_v2_genomic_2500GB.bed  GCF_000233375.1_ICSASG_v2_genomic_2500GE.bed GCF_000233375.1_ICSASG_v2_genomic_GB.bed ; do  bedtools intersect -a dmrs.bed -b $feature -wb | awk -v i="$feature" '{print $16, i}'| tr ";" "\t" | awk '{print $2, $NF}'| cut -c 15- >> genes_X_DMR.txt ; done
for feature in GCF_000233375.1_ICSASG_v2_genomic_250promoter.bed GCF_000233375.1_ICSASG_v2_genomic_750promoter.bed GCF_000233375.1_ICSASG_v2_genomic_5000promoter.bed GCF_000233375.1_ICSASG_v2_genomic_2500GB.bed  GCF_000233375.1_ICSASG_v2_genomic_2500GE.bed GCF_000233375.1_ICSASG_v2_genomic_GB.bed ; do  bedtools intersect -a dmrs.bed -b $feature -wb | awk -v i="$feature" '{print $10, i}' |cut -c 6- >> DMR_X_genes.txt ; done

## Run in R
ncbi<-read_tsv("~/Downloads/gene_result.txt") %>% dplyr::select(GeneID, Symbol, description, Aliases,Status)
dmrs<-read_tsv("~/Documents/Crappyfish/DMR_X_genes.txt", col_names = c("seqnames", "start", "end",".X1","X2","feature"))
genes<-read_delim("~/Documents/Crappyfish/genes_X_DMR.txt", delim=" ", col_names = c("gene", "feature"))

dmr_genes<-tibble(dmrs, "gene"=genes$gene)
table<-as_tibble(mydmr) %>% 
  rowid_to_column(.,"DMR") %>% 
  full_join(dmr_genes,., by=c("seqnames"="seqnames","end"="end")) %>% 
  dplyr::select(DMR,"chromosome"=seqnames,"start"=start.y,end, mean.meth.diff,num.CpGs,feature,gene) %>% 
  arrange(chromosome,start) %>% 
  mutate(feature = replace(feature, feature == ". GCF_000233375.1_ICSASG_v2_genomic_5000promoter.bed", "Distal promoter"))  %>% 
  mutate(feature = replace(feature, feature == ". GCF_000233375.1_ICSASG_v2_genomic_750promoter.bed", "Medial promoter"))  %>% 
  mutate(feature = replace(feature, feature == ". GCF_000233375.1_ICSASG_v2_genomic_250promoter.bed", "Proximal promoter"))  %>% 
  mutate(feature = replace(feature, feature == ". GCF_000233375.1_ICSASG_v2_genomic_2500GE.bed", "Gene end")) %>% 
  mutate(feature = replace(feature, feature == ". GCF_000233375.1_ICSASG_v2_genomic_2500GB.bed", "Gene start")) %>%
  mutate(feature = replace(feature, feature == ". GCF_000233375.1_ICSASG_v2_genomic_GB.bed", "Gene body")) %>%
  replace_na(list(feature="Intergenic")) %>%
  unique()  

tableS2<-left_join(table,ncbi,by=c("gene"="GeneID"))
view(tableS2)
