library("DESeq2")
library("tidyverse")

setwd("/home/benjamin/Documents/Brassicas_repo")

#oull relevant metadata
data.meta = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info.csv") %>%
  mutate(Domesticated = ifelse(Wild..Domesticated=="Wild","Wild","Cultivated")) %>%
  mutate(Environment = ifelse(Environment=="wheat competition","wheat","control"),
         Parental.effects.status = ifelse(Parental.effects.status == "'\"standardised\"","standardised","unstandardised")) %>%
  select(c("RNAseq.sample.name","Species","Parental.effects.status","Environment","Domesticated")) %>%
  'colnames<-'(c("sample","species","parental.effects","treatment","domesticated")) %>%
  mutate_all(as.factor)

#subset by genus
metadata.brass = subset(data.meta, substr(species,1,3)=="Bra")
metadata.raph = subset(data.meta, substr(species,1,3)=="Rap")

#pull rnaseq data
raph.gene.counts = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/raph.gene.counts.csv", row.names = 1)



#filter by expression
