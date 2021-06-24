library("DESeq2")
library("tidyverse")
library("WGCNA")

setwd("/home/benjamin/Documents/Brassicas_repo")

#import functions
source("Functions/expression_heatmap.R")
source("Functions/topGO_wrapper.R")

#pull relevant metadata
data.meta = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info.csv") %>%
  mutate(Domesticated = ifelse(Wild..Domesticated=="Wild","Wild","Cultivated")) %>%
  mutate(Environment = ifelse(Environment=="wheat competition","Wheat","Control"),
         Parental.effects.status = ifelse(Parental.effects.status == "'\"standardised\"","Standardised","Unstandardised")) %>%
  dplyr::select(c("RNAseq.sample.name","Species","Parental.effects.status","Environment","Domesticated")) %>%
  #dplyr::select(c("RNAseq.sample.name","Species","Parental.effects.status","Environment","Wild..Domesticated")) %>%
  'colnames<-'(c("sample","species","parental.effects","treatment","domesticated")) %>%
  mutate_all(as.factor)

data.meta = mutate(data.meta, 
                   treatment = fct_relevel(treatment, c("Control","Wheat")),
                   domesticated = fct_relevel(domesticated, c("Wild","Cultivated")))


#subset by genus
metadata.brass = subset(data.meta, substr(species,1,3)=="Bra")
metadata.raph = subset(data.meta, substr(species,1,3)=="Rap")

#pull rnaseq data
brass.gene.counts = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/brapa.gene.counts.csv", row.names = 1)

#filter by expression
brass.gene.counts.clean = brass.gene.counts[(which(rowMeans(brass.gene.counts)>=1)),]
paste0(nrow(brass.gene.counts)-nrow(brass.gene.counts.clean),"/",nrow(brass.gene.counts)," Brassica genes filtered due to very low expression.")

#an additional filter- remove genes that have too many 0 counts
paste0("Removing an additional ",nrow(brass.gene.counts.clean[rowSums(brass.gene.counts.clean >= 5) < 3,])," genes with many low counts.")
brass.gene.counts.clean = brass.gene.counts.clean[rowSums(brass.gene.counts.clean >= 5) >= 3,]

#check that metadata and count matrices conform
table(metadata.brass$sample %in% colnames(brass.gene.counts.clean))
brass.gene.counts.clean = brass.gene.counts.clean[,as.character(metadata.brass$sample)]

#some basic QC: heatmaps and pca
heatMap = expression.heatmap(countdata = sample_n(brass.gene.counts.clean,100),
                             data.phenotype = metadata.brass,
                             labels = c("parental.effects",
                                        "domesticated"
                                        #,
                             #            "treatment",
                             #            "domesticated"
                             ),
                             pass_on = F,
                             ID_var = "sample")

# log-transform with a pseudocount
pca.counts = log2(brass.gene.counts.clean+1)
#create pca object
data.pca = prcomp(t(pca.counts))
#extract PC data
percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
pca.out = list(values = data.frame(data.pca$x),
               percent.var = percent.var)
#connect to phenotypic data
ggpcadata = pca.out$values %>%
  rownames_to_column(var = "sample") %>%
left_join(metadata.brass,
            by = "sample")
#plot
ggplot(ggpcadata, aes(x = PC1, y = PC2, shape = domesticated, color = parental.effects, label = sample)) +
  geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
  #geom_text(vjust = -1) +
  xlab(paste0("PC",1,": ",signif(pca.out$percent.var[1]*100, 3),"%")) +
  ylab(paste0("PC",2,": ",signif(pca.out$percent.var[2]*100, 3),"%")) +
  theme_bw() +
  # scale_color_manual(name = "Treatment",
  #                    values = brewer.pal(7, "Paired")) +
  scale_shape_manual(name = "Treatment",
                     values = c(8,15:20)) +
  theme(panel.grid = element_line(color = "grey95"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(face = "bold", size =12))


#next step is to check for genes with parental effects and exclude these
dds.parental = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean,
                                      colData = metadata.brass,
                                      design = as.formula(~parental.effects))
dds.parental.deg = DESeq(dds.parental, fitType = "local", betaPrior = FALSE)
parental.degs = results(dds.parental.deg)
#a high proportion of genes have parental effects
print(paste0("Number of genes with parental effects at p<0.1: ",length(which(parental.degs$padj<0.1)),"/",nrow(parental.degs)))
parental.degs.ids = rownames(subset(parental.degs, padj<=0.1))
#brass.gene.counts.clean = brass.gene.counts.clean[rownames(subset(parental.degs, padj>=0.1)),]
#also check for parental deg GO terms 
GOscores.parental = as.numeric(row.names(brass.gene.counts.clean)%in%parental.degs.ids) %>% 'names<-'(row.names(brass.gene.counts.clean))
parentGO = topGO_wrapper(geneScores = GOscores.parental,
                         geneScoresDE = F,
                         geneScoresDirection = NA,
                         GOmapping = GOmapping.brass,
                         algorithm = "weight01",
                         statistic = "fisher",
                         nodeSize = 10,
                         discretisedDE = F,
                         p = 0.05)
parentGO$consolidated_result
#47 GO terms, many biosynthetic and developmental 

#now run proper model
dds.gene = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean,
                                  colData = metadata.brass,
                                  design = as.formula(~treatment+domesticated+treatment*domesticated))

# Run the default analysis for DESeq2 and generate results table. NA p-values are generated by 0 counts and outliers calculated by Cook's distance.
dds.gene.deg = DESeq(dds.gene, fitType = "local", betaPrior = FALSE)
# Check for outliers: none are apparent
print("Check for gene expression outliers")
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds.gene.deg)[["cooks"]]), range=0, las=2)

#Also plot the dispersion extimates to make sure they look fine
plotDispEsts(dds.gene.deg)

#check effectof a parametric rather than local fit
parametric.gene.deg = DESeq(dds.gene, fitType = "parametric", betaPrior = FALSE)
plotDispEsts(parametric.gene.deg) 
# parametric seems to follow the data a little better so let's go with that
dds.gene.deg = parametric.gene.deg

resultsNames(dds.gene.deg)

#wheat vs control
degs.brass.treatment = results(dds.gene.deg, 
                              name="treatment_Wheat_vs_Control",     
                              alpha = 0.05,
                              lfcThreshold = log2(1))
summary(degs.brass.treatment) #15 DEGs
degs.brass.treatment.ids = rownames(subset(degs.brass.treatment, padj<=0.05))
#1 of these degs is shared with parental degs
table(degs.brass.treatment.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
#GOscores.brass.treatment = as.numeric(row.names(brass.gene.counts.clean)%in%degs.brass.treatment.ids) %>% 'names<-'(row.names(brass.gene.counts.clean))
brass.treatment.GO.up = topGO_wrapper(geneScores = degs.brass.treatment,
                                  geneScoresDE = T,
                                  geneScoresDirection = "Up",
                                  GOmapping = GOmapping.brass,
                                  algorithm = "weight01",
                                  statistic = "fisher",
                                  nodeSize = 10,
                                  discretisedDE = T,
                                  p = 0.05)
brass.treatment.GO.up$consolidated_result
#25 GO terms: mostly metabolic but also a couple defensive against fungus e.g. "defense response to fungus" and "response to chitin"
#"also response to cold" and "response to water deprivation"
brass.treatment.GO.down = topGO_wrapper(geneScores = degs.brass.treatment,
                                      geneScoresDE = T,
                                      geneScoresDirection = "Down",
                                      GOmapping = GOmapping.brass,
                                      algorithm = "weight01",
                                      statistic = "fisher",
                                      nodeSize = 10,
                                      discretisedDE = T,
                                      p = 0.05)
brass.treatment.GO.down$consolidated_result


#cultivated vs wild
degs.brass.cultivated = results(dds.gene.deg, 
                               name="domesticated_Cultivated_vs_Wild",     
                               alpha = 0.05,
                               lfcThreshold = log2(1))
summary(degs.brass.cultivated) #~1500 degs
degs.brass.cultivated.ids = rownames(subset(degs.brass.cultivated, padj<=0.05))
#12 of these degs is shared with parental degs
table(degs.brass.cultivated.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
#GOscores.brass.cultivated = as.numeric(row.names(brass.gene.counts.clean)%in%degs.brass.cultivated.ids) %>% 'names<-'(row.names(brass.gene.counts.clean))
brass.cultivated.GO.up = topGO_wrapper(geneScores = degs.brass.cultivated,
                                   geneScoresDE = T,
                                   geneScoresDirection = "Up",
                                   GOmapping = GOmapping.brass,
                                   algorithm = "weight01",
                                   statistic = "fisher",
                                   nodeSize = 10,
                                   discretisedDE = T,
                                   p = 0.05)
brass.cultivated.GO.up$consolidated_result
#34 GO terms: lots of developmental, cell growth and death, defensive against fungus, 

brass.cultivated.GO.down = topGO_wrapper(geneScores = degs.brass.cultivated,
                                       geneScoresDE = T,
                                       geneScoresDirection = "Down",
                                       GOmapping = GOmapping.brass,
                                       algorithm = "weight01",
                                       statistic = "fisher",
                                       nodeSize = 10,
                                       discretisedDE = T,
                                       p = 0.05)
brass.cultivated.GO.down$consolidated_result
#52 terms, many related to cell cycle and DNA replication, chromatin silencing, gene silcing by mirna


#interaction
degs.brass.interaction = results(dds.gene.deg, 
                                name="treatmentWheat.domesticatedCultivated",     
                                alpha = 0.05,
                                lfcThreshold = log2(1))
summary(degs.brass.interaction) #46 degs
degs.brass.interaction.ids = rownames(subset(degs.brass.interaction, padj<=0.05))
#1 of these degs is shared with parental degs
table(degs.brass.interaction.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
#GOscores.brass.interaction = as.numeric(row.names(brass.gene.counts.clean)%in%degs.brass.interaction.ids) %>% 'names<-'(row.names(brass.gene.counts.clean))
brass.interaction.GO = topGO_wrapper(geneScores = degs.brass.interaction,
                                    geneScoresDE = T,
                                    geneScoresDirection = NA,
                                    GOmapping = GOmapping.brass,
                                    algorithm = "weight01",
                                    statistic = "fisher",
                                    nodeSize = 10,
                                    discretisedDE = T,
                                    p = 0.05)
brass.interaction.GO$consolidated_result
#26 GO terms: many epigenetic "histone binding" and "histone methylation", "chromatin binding", "regulation of histone modification"
#several response to sugar (fructose, glucose, sucrose), many developmental
#for the interaction terms, we also plot the output to understand what exactly is going on
#get the 12 terms with lowest padj
degs.brass.interaction.signif = subset(degs.brass.interaction, padj<0.05)
interestgenes = row.names(degs.brass.interaction.signif)[order(degs.brass.interaction.signif$log2FoldChange,decreasing = T)[1:12]]
#for each gene, extract the counts for plotting and label by gene
for(i in 1:length(interestgenes)){
  #if first element, instantiate frame, otherwise rbind to frame
  if(i==1) {
    brass.intplotdata = plotCounts(dds.gene.deg, gene=interestgenes[i], intgroup=c("domesticated","treatment"), returnData = T)
    brass.intplotdata$gene = interestgenes[i]
  } else {
    addrows = data.frame(plotCounts(dds.gene.deg, gene=interestgenes[i], intgroup=c("domesticated","treatment"), returnData = T), 
                         gene = interestgenes[i])
    brass.intplotdata = rbind(brass.intplotdata, addrows)
  }
}
#plot with ggplot facet wrap
brass.intplot= ggplot(brass.intplotdata, aes(x = treatment, y = count)) +
  stat_summary(aes(group = domesticated), fun.y = mean, geom = "path") +
  stat_summary(aes(color = domesticated), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(color = domesticated), fun.y = mean, geom = "point", size = 4) +
  geom_point(aes(color = domesticated), size = 2) +
  scale_color_manual(values = brewer.pal(6,"Set3")[c(5,6)]) +
  facet_wrap(~gene, scales = "free")

ggsave(brass.intplot, 
       filename = "brass_interaction_deg_plots.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  40, height = 25, units = "cm")