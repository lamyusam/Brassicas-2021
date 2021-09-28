setwd("/home/benjamin/Documents/Brassicas_repo")

# get libraries
basic_libraries = c("DESeq2",
                    "tidyverse",
                    "WGCNA",
                    "RColorBrewer",
                    "topGO",
                    "rrvgo",
                    "GOSemSim",
                    "org.Dm.eg.db")
for(lib in basic_libraries){
  if(require(package = lib, character.only = TRUE)){
    print("Successful")
  }else{print("Installing")
    install.packages(lib, Ncpus = 6)
    library(lib, character.only = TRUE)}
}
#prevent other packages overriding dplyr's select
select = dplyr::select
#import custom functions
source("Functions/expression_heatmap.R")
source("Functions/topGO_wrapper.R")
source("Functions/gene_overlap_test.R")
source("Functions/GO_treeplots.R")

## Data import, pre-filtering and QC
#pull relevant metadata
data.meta.gen4 = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info_gen4.csv") %>%
  #rename variable levels for consistency
  mutate(Domesticated = ifelse(Wild..Domesticated=="Wild","Wild","Cultivated")) %>%
  mutate(Environment = ifelse(Environment=="wheat competition","Wheat","Control"),
         Parental.effects.status = ifelse(Parental.effects.status == "'\"standardised\"","Standardised","Unstandardised")) %>%
  #drop unused variables
  dplyr::select(c("RNAseq.sample.name","Species","Label","Parental.effects.status","Environment","Domesticated","Generation")) %>%
  #rename variables for consistency
  'colnames<-'(c("sample","species","label","parental.effects","treatment","domesticated","generation")) %>%
  #convert all variables to factor
  mutate_all(as.factor) %>%
  #relevel variables of interest so that the 'natural' state (wild+wheat) is always the base level
  mutate(treatment = fct_relevel(treatment, c("Wheat","Control")),
         domesticated = fct_relevel(domesticated, c("Wild","Cultivated")))

#subset by genus
metadata.brass.gen4 = subset(data.meta.gen4, substr(species,1,3)=="Bra")
metadata.raph.gen4 = subset(data.meta.gen4, substr(species,1,3)=="Rap")

#import rnaseq data
brass.gene.counts.gen4 = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/brapa.gene.counts.gen4.csv", row.names = 1)
raph.gene.counts.gen4 = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/raph.gene.counts.gen4.csv", row.names = 1)
#import GO mapping
load("Data/GO/brass_GOmapping.Rdata")
load("Data/GO/raph_GOmapping.Rdata")


#filter by expression
brass.gene.counts.gen4.clean = brass.gene.counts.gen4[(which(rowMeans(brass.gene.counts.gen4)>=1)),]
print(paste0(nrow(brass.gene.counts.gen4)-nrow(brass.gene.counts.gen4.clean),"/",nrow(brass.gene.counts.gen4)," Brassica genes filtered due to very low expression."))
#an additional filter- remove genes that have too many 0 counts
print(paste0("Removing an additional ",nrow(brass.gene.counts.gen4.clean[rowSums(brass.gene.counts.gen4.clean >= 5) < 3,])," Brassica genes with many low counts."))
brass.gene.counts.gen4.clean = brass.gene.counts.gen4.clean[rowSums(brass.gene.counts.gen4.clean >= 5) >= 3,]
#ensure that metadata and count matrices conform
metadata.brass.gen4 = metadata.brass.gen4[which((metadata.brass.gen4$sample %in% colnames(brass.gene.counts.gen4.clean))),]
brass.gene.counts.gen4.clean = brass.gene.counts.gen4.clean[,as.character(metadata.brass.gen4$sample)]
# log-transform with a pseudocount for PCA
pca.counts = log2(brass.gene.counts.gen4.clean+1)
#create pca object
data.pca = prcomp(t(pca.counts))
#extract PC data
percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
pca.out = list(values = data.frame(data.pca$x),
               percent.var = percent.var)
#connect to phenotypic data
ggpcadata = pca.out$values %>%
  rownames_to_column(var = "sample") %>%
  left_join(metadata.brass.gen4,
            by = "sample")
#plot but with species labels
brass.pca = ggplot(ggpcadata, aes(x = PC1, y = PC2, shape = domesticated, label = sample)) +
  geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
  geom_text(vjust = -1) +
  xlab(paste0("PC",1,": ",signif(pca.out$percent.var[1]*100, 3),"%")) +
  ylab(paste0("PC",2,": ",signif(pca.out$percent.var[2]*100, 3),"%")) +
  theme_bw() +
  # scale_color_manual(name = "Treatment",
  #                    values = brewer.pal(7, "Paired")) +
  scale_shape_manual(name = "Domestication",
                     values = c(8,15:20)) +
  theme(panel.grid = element_line(color = "grey95"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(face = "bold", size =12))
#save plot
ggsave(brass.pca, 
       filename = "brass_clean_pca_gen4.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  25, height = 15, units = "cm")
#plot PCA for markdown
knitr::include_graphics("Analysis/RNAseq/Images/brass_clean_pca_gen4.png")

#from this we can see odd clustering of b rapas
checkframe = read.csv("Data/RNAseq/RNASeq_sample_info_gen4.csv")
#the brapa samples clearly divide into two very distinct groups along the first PC
outgroup = subset(ggpcadata, PC1>(100))$sample
outcheck= subset(checkframe, RNAseq.sample.name %in% outgroup)
#again, the outleirs are all tricoloris.... we'll keep this in for now but it's a pain
tricoloris.outsamples = subset(checkframe, RNAseq.sample.name %in% outgroup)$RNAseq.sample.name


## Now repeat QC and pre-filtering for Raphanus samples:
  
#filter by expression, removing genes with <1 count/sample
raph.gene.counts.gen4.clean = raph.gene.counts.gen4[(which(rowMeans(raph.gene.counts.gen4)>=1)),]
paste0(nrow(raph.gene.counts.gen4)-nrow(raph.gene.counts.gen4.clean),"/",nrow(raph.gene.counts.gen4)," Raphanus genes filtered due to very low expression.")
#an additional filter- remove genes that only have expression in a couple of samples
paste0("Removing an additional ",nrow(raph.gene.counts.gen4.clean[rowSums(raph.gene.counts.gen4.clean >= 5) < 3,])," genes with many low counts.")
raph.gene.counts.gen4.clean = raph.gene.counts.gen4.clean[rowSums(raph.gene.counts.gen4.clean >= 5) >= 3,]
#check that metadata and count matrices conform
#table(metadata.raph.gen4$sample %in% colnames(raph.gene.counts.gen4.clean))
raph.gene.counts.gen4.clean = raph.gene.counts.gen4.clean[,as.character(metadata.raph.gen4$sample)]

# log-transform with a pseudocount for PCA
pca.counts = log2(raph.gene.counts.gen4.clean+1)
#create pca object
data.pca = prcomp(t(pca.counts))
#extract PC data
percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
pca.out = list(values = data.frame(data.pca$x),
               percent.var = percent.var)
#connect to phenotypic data
ggpcadata = pca.out$values %>%
  rownames_to_column(var = "sample") %>%
  left_join(metadata.raph.gen4,
            by = "sample")
#plot
raph.clean.pca.gen4 = ggplot(ggpcadata, aes(x = PC1, y = PC2, color = species, shape = domesticated, label = sample)) +
  geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
  #geom_text(vjust = -1) +
  xlab(paste0("PC",1,": ",signif(pca.out$percent.var[1]*100, 3),"%")) +
  ylab(paste0("PC",2,": ",signif(pca.out$percent.var[2]*100, 3),"%")) +
  theme_bw() +
  scale_shape_manual(name = "Domestication",
                     values = c(8,15:20)) +
  theme(panel.grid = element_line(color = "grey95"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(face = "bold", size =12))
# #save plot
ggsave(raph.clean.pca.gen4,
       filename = "raph_clean_pca_gen4.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  25, height = 15, units = "cm")
#We see unusual clustering in the PCA for Raphanus, as with Brassica. The 'outlying' samples on PC2 are represent all Asian sativus var. caudatus, while the U.S. caudatus samples cluster with the other wilds. As with the unusual B. tricoloris samples, we'll retain these for now but make a note in case we wish to exclude them at a later time. 
#check samples that look unusual
checkframe = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info_gen4.csv") 
subset(checkframe , Label %in% subset(ggpcadata,PC2<(-100))$label)
outgroup = subset(ggpcadata,PC2<(-100))$label
caudatus.outsamples = subset(checkframe, Label %in% outgroup)$RNAseq.sample.name
knitr::include_graphics("Analysis/RNAseq/Images/raph_clean_pca.png")

