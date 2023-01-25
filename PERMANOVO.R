library(vegan)

load("{your_path}/MethylKit_19samples.RData")

sample.id <- c("D1","D2","D3","D4","D5","D6","D8","D9","D10", "D11","D12","D13","D14","D15","D16","D17","D18","D19","D20")
treatment_n <- c(1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0)
treatment <- gsub(1,"Sick",treatment_n) %>% gsub("0","Healthy",.)
design.matrix <- bind_cols(ID=sample.id, treatment=treatment)


PercM <- percMethylation(meth.filtered) %>% t(.) %>% as_tibble(.)

# Perform permanova test using adonis function
a<-adonis2(PercM ~ treatment, method = "euclidean")
a
