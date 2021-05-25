library(tidyverse)
library(ggbeeswarm)
#library(ggforce)

setwd("/home/benjamin/Documents/Brassicas_repo/")

#get phenotypic data for later
phenodata.gen2 = read.delim("Data/Phenotypic_data/morphology_data_gen2.txt", sep = ",", header = T, row.names = NULL)

#get alignment data for brassica
alignrates.brass = read.delim("Analysis/RNAseq/QC/alignment_rates/alignmentdata_brassica.tsv", header = F, sep = " ") %>%
  'colnames<-'(c("sample","reads_total","reads_aligned")) %>% 
  mutate(sample = as.numeric(stringr::str_match(sample,"[0-9]+")))
#correct labels (necessary because I messed up the nextflow sample sheet)
nextflowsheet.brass = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/nextflow_samplesheet_brass.csv") %>%
  mutate(sampleID = stringr::str_match(fastq_1,"A[0-9]+")[,1])
alignrates.brass$sample = nextflowsheet.brass$sampleID[match(alignrates.brass$sample,nextflowsheet.brass$replicate)]

#get alignment data for raphanus
alignrates.raph = read.delim("Analysis/RNAseq/QC/alignment_rates/alignmentdata_raphanus.tsv", header = F, sep = " ") %>%
  'colnames<-'(c("sample","reads_total","reads_aligned")) %>% 
  mutate(sample = as.numeric(stringr::str_match(sample,"[0-9]+")))
#correct labels (necessary because I missed up w the nextflow sample sheet)
nextflowsheet.raph = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/nextflow_samplesheet_raph.csv") %>%
  mutate(sampleID = stringr::str_match(fastq_1,"A[0-9]+")[,1])
alignrates.raph$sample = nextflowsheet.raph$sampleID[match(alignrates.raph$sample,nextflowsheet.raph$replicate)]

#consolidate alignment data
alignrates = rbind(alignrates.brass,alignrates.raph) %>%
  mutate(alignrate = (reads_aligned/reads_total)*100)
  
#get rnaseq metadata and attach alignment data
rna.metadata.gen2 = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info.csv")[,-1] %>%
  left_join(alignrates[,c("sample","alignrate")],by=c("RNAseq.sample.name" = "sample")) %>%
  mutate(SpeciesShort = Species)

#shorten names for plotting
rna.metadata.gen2$SpeciesShort = gsub("Brassica","B.",rna.metadata.gen2$SpeciesShort)
rna.metadata.gen2$SpeciesShort = gsub("Raphanus","R.",rna.metadata.gen2$SpeciesShort)

#plot by species
ggplot(rna.metadata.gen2, aes(x = substr(Species,1,4), y = alignrate)) +
  #geom_violin() +
  geom_beeswarm() +
  xlab("Species") +
  ylab("Alignment rate (%)") +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() + 
  theme(aspect.ratio = 2,
        panel.grid = element_line(size = 0.2, colour = "gray80"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size =15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=15, face = "bold")) 

#ggplot(subset(rna.metadata.gen2,substr(Species,1,4)=="Bras"), aes(x = Wild..Domesticated, y = alignrate)) +
brassica.byspecies = ggplot(subset(rna.metadata.gen2,substr(Species,1,4)=="Bras"), aes(x = SpeciesShort, y = alignrate)) +
  #geom_violin(size = 0.3, fill = "gray95") +
  #geom_beeswarm() +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  xlab("Species") +
  ylab("Alignment rate (%)") +
  scale_y_continuous(limits = c(60,100)) +
  theme_bw() + 
  theme(aspect.ratio = 0.8,
        panel.grid = element_line(size = 0.2, colour = "gray80"),
        axis.text.x = element_text(size = 12, angle = 0),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size =15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=15, face = "bold")) 

#save
ggsave(brassica.byspecies, filename = "brassica_alignment_byspecies.png", device = "png", 
       path ="Analysis/RNAseq/QC/alignment_rates", 
       width = 30, height = 18, units = "cm")

brassica.byspecies = ggplot(subset(rna.metadata.gen2,substr(Species,1,4)=="Raph"), aes(x = SpeciesShort, y = alignrate)) +
  #geom_violin(size = 0.3, fill = "gray95") +
  #geom_beeswarm() +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  xlab("Species") +
  ylab("Alignment rate (%)") +
  scale_y_continuous(limits = c(60,100)) +
  theme_bw() + 
  theme(aspect.ratio = 0.8,
        panel.grid = element_line(size = 0.2, colour = "gray80"),
        axis.text.x = element_text(size = 10, angle = 0),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size =15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=15, face = "bold")) 

#save
ggsave(brassica.byspecies, filename = "raphanus_alignment_byspecies.png", device = "png", 
       path ="Analysis/RNAseq/QC/alignment_rates", 
       width = 30, height = 18, units = "cm")

