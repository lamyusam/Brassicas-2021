# get libraries
basic_libraries <- c("DESeq2",
                     "tidyverse",
                     "WGCNA",
                     "RColorBrewer")

for (lib in basic_libraries) {
  if (require(package = lib, character.only = TRUE)) {
    print("Successful")
  } else {
    print("Installing")
    install.packages(lib)
    library(lib, character.only = TRUE )
  }
}

#set working directory
setwd("/home/benjamin/Documents/Brassicas_repo")

#import custom functions
source("Functions/expression_heatmap.R")
source("Functions/topGO_wrapper.R")

#import GO mapping
load("Data/GO/raph_GOmapping.Rdata")

#pull relevant metadata and rename/relevel variables
data.meta = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info.csv") %>%
  mutate(Domesticated = ifelse(Wild..Domesticated=="Wild","Wild","Cultivated")) %>%
  mutate(Environment = ifelse(Environment=="wheat competition","Wheat","Control"),
         Parental.effects.status = ifelse(Parental.effects.status == "'\"standardised\"","Standardised","Unstandardised")) %>%
  select(c("RNAseq.sample.name","Species","Parental.effects.status","Environment","Domesticated")) %>%
  'colnames<-'(c("sample","species","parental.effects","treatment","domesticated")) %>%
  mutate_all(as.factor)

data.meta = mutate(data.meta, 
                   treatment = fct_relevel(treatment, c("Control","Wheat")),
                   domesticated = fct_relevel(domesticated, c("Wild","Cultivated")))


#subset by genus
metadata.brass = subset(data.meta, substr(species,1,3)=="Bra")
metadata.raph = subset(data.meta, substr(species,1,3)=="Rap")

#save for mark
write.csv(metadata.raph, file = "Data/RNAseq/metadata_raph.csv", row.names = F)

#pull rnaseq data
raph.gene.counts = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/raph.gene.counts.csv", row.names = 1)

#filter by expression, removing genes with <1 count/sample
raph.gene.counts.clean = raph.gene.counts[(which(rowMeans(raph.gene.counts)>=1)),]
paste0(nrow(raph.gene.counts)-nrow(raph.gene.counts.clean),"/",nrow(raph.gene.counts)," Raphanus genes filtered due to very low expression.")

#an additional filter- remove genes that only have expression in a couple of samples
paste0("Removing an additional ",nrow(raph.gene.counts.clean[rowSums(raph.gene.counts.clean >= 5) < 3,])," genes with many low counts.")
raph.gene.counts.clean = raph.gene.counts.clean[rowSums(raph.gene.counts.clean >= 5) >= 3,]

#check that metadata and count matrices conform
table(metadata.raph$sample %in% colnames(raph.gene.counts.clean))
raph.gene.counts.clean = raph.gene.counts.clean[,as.character(metadata.raph$sample)]

#some basic QC: heatmaps and pca
heatMap = expression.heatmap(countdata = sample_n(raph.gene.counts.clean,5000),
                             data.phenotype = metadata.raph,
                             labels = c("parental.effects",
                                        "treatment",
                                        "domesticated"),
                             pass_on = F,
                             ID_var = "sample")

# log-transform with a pseudocount
pca.counts = log2(raph.gene.counts.clean+1)
#create pca object
data.pca = prcomp(t(pca.counts))
#extract PC data
percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
pca.out = list(values = data.frame(data.pca$x),
               percent.var = percent.var)
#connect to phenotypic data
ggpcadata = pca.out$values %>%
  rownames_to_column(var = "sample") %>%
  left_join(metadata.raph,
            by = "sample")
#plot
ggplot(ggpcadata, aes(x = PC1, y = PC2, color = domesticated, shape = treatment, label = sample)) +
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

#next step is to check for genes with parental effects
dds.parental = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean,
                                      colData = metadata.raph,
                                      design = as.formula(~parental.effects))
dds.parental.deg = DESeq(dds.parental, fitType = "parametric", betaPrior = FALSE)
parental.degs = results(dds.parental.deg)
#very few genes have even marginal evidence of parental effects
print(paste0("Number of genes with parental effects at p<0.1: ",length(which(parental.degs$padj<0.1)),"/",nrow(parental.degs)))
parental.degs.ids = rownames(subset(parental.degs, padj<=0.1))
#uncomment next line to simply exclude these genes from the analysis
#raph.gene.counts.clean = raph.gene.counts.clean[rownames(subset(parental.degs, padj>=0.1)),]
#also check for parental deg GO terms 
parentGO = topGO_wrapper(geneScores = parental.degs,
                         geneScoresDE = T,
                         geneScoresDirection = NA,
                         GOmapping = GOmapping.raph,
                         algorithm = "weight01",
                         statistic = "fisher",
                         nodeSize = 10,
                         discretisedDE = T,
                         p = 0.05)
parentGO$consolidated_result
#48 GO terms- a surprisingly high number given the relatively small number of genes 

#subset to remove all but raphanistrum sativus and its wild ancestor
metadata.raph.subset = subset(metadata.raph, species %in% c("Raphanus raphanistrum","Raphanus sativus"))
raph.gene.counts.clean.subset = raph.gene.counts.clean[,as.character(metadata.raph.subset$sample)]
#check that metadata and count matrices conform
table(metadata.raph.subset$sample %in% colnames(raph.gene.counts.clean.subset))
raph.gene.counts.clean.subset = raph.gene.counts.clean.subset[,as.character(metadata.raph.subset$sample)]

#now run model with our variables of interest
dds.gene = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean.subset,
                                  colData = metadata.raph.subset,
                                  design = as.formula(~treatment+domesticated+treatment*domesticated))

# Run the default analysis for DESeq2 and generate results table. 
# NA p-values are generated by 0 counts and outliers calculated by Cook's distance.
dds.gene.deg = DESeq(dds.gene, fitType = "local", betaPrior = FALSE)
# Check for outliers: none are apparent
print("Check for gene expression outliers")
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds.gene.deg)[["cooks"]]), range=0, las=2)

#Also plot the dispersion extimates to make sure they look fine
plotDispEsts(dds.gene.deg)

#check effect of a parametric rather than local fit
parametric.gene.deg = DESeq(dds.gene, fitType = "parametric", betaPrior = FALSE)
plotDispEsts(parametric.gene.deg) 

# parametric seems to follow the data better so we'll use that going forward
dds.gene.deg = parametric.gene.deg

#results: wheat vs control
degs.raph.treatment = results(dds.gene.deg, 
                              name="treatment_Wheat_vs_Control",     
                              alpha = 0.05,
                              lfcThreshold = log2(1))
summary(degs.raph.treatment) #126 DEGs
degs.raph.treatment.ids = rownames(subset(degs.raph.treatment, padj<=0.05))
#1 of these degs is shared with parental degs
table(degs.raph.treatment.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
raph.treatment.GO.up = topGO_wrapper(geneScores = degs.raph.treatment,
                         geneScoresDE = T,
                         geneScoresDirection = "Up",
                         GOmapping = GOmapping.raph,
                         algorithm = "weight01",
                         statistic = "fisher",
                         nodeSize = 10,
                         discretisedDE = T,
                         p = 0.05)
raph.treatment.GO.up$consolidated_result
#21 GO terms up, inc. histone binding, chromatin binding, 'response to chitin', response to salt stress, 
#response to sucrose, response to glucose, plant ovule development, regulation of flower development, telomere maintenance
#terms down:
raph.treatment.GO.down = topGO_wrapper(geneScores = degs.raph.treatment,
                                     geneScoresDE = T,
                                     geneScoresDirection = "Down",
                                     GOmapping = GOmapping.raph,
                                     algorithm = "weight01",
                                     statistic = "fisher",
                                     nodeSize = 10,
                                     discretisedDE = T,
                                     p = 0.05)
raph.treatment.GO.down$consolidated_result
#19 Go terms down, inc. cold acclimation, response to cold, response to chitin, response to water deprivation, response to ozone 

#results: cultivated vs wild
degs.raph.cultivated = results(dds.gene.deg, 
                               name="domesticated_Cultivated_vs_Wild",     
                               alpha = 0.05,
                               lfcThreshold = log2(1))
summary(degs.raph.cultivated) #~1500 degs
degs.raph.cultivated.ids = rownames(subset(degs.raph.cultivated, padj<=0.05))
#12 of these degs is shared with parental degs
table(degs.raph.cultivated.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
#GOscores.raph.cultivated = as.numeric(row.names(raph.gene.counts.clean.subset)%in%degs.raph.cultivated.ids) %>% 'names<-'(row.names(raph.gene.counts.clean.subset))
raph.cultivated.GO.up = topGO_wrapper(geneScores = degs.raph.cultivated,
                                  geneScoresDE = T,
                                  geneScoresDirection = "Up",
                                  GOmapping = GOmapping.raph,
                                  algorithm = "weight01",
                                  statistic = "fisher",
                                  nodeSize = 10,
                                  discretisedDE = T,
                                  p = 0.05)
raph.cultivated.GO.up$consolidated_result
#33 GO terms up: mostly metabolic and developmental but also response to cold, chromatin assembly, aging, 
#leaf  senescence, regulation of chromatin organization
raph.cultivated.GO.down = topGO_wrapper(geneScores = degs.raph.cultivated,
                                      geneScoresDE = T,
                                      geneScoresDirection = "Down",
                                      GOmapping = GOmapping.raph,
                                      algorithm = "weight01",
                                      statistic = "fisher",
                                      nodeSize = 10,
                                      discretisedDE = T,
                                      p = 0.05)
raph.cultivated.GO.down$consolidated_result
#57 down: heavy on the chloroplasts, 'choroplast','chloroplast envelope','photorespriartion','response to chitin', 'response to chitin',
#'response to woudning', lots and lots of biosynthetic processes

#results: interaction
degs.raph.interaction = results(dds.gene.deg, 
                                name="treatmentWheat.domesticatedCultivated",     
                                alpha = 0.05,
                                lfcThreshold = log2(1))
summary(degs.raph.interaction) #96 degs
degs.raph.interaction.ids = rownames(subset(degs.raph.interaction, padj<=0.05))
#1 of these degs is shared with parental degs
table(degs.raph.interaction.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
GOscores.raph.interaction = as.numeric(row.names(raph.gene.counts.clean.subset)%in%degs.raph.interaction.ids) %>% 'names<-'(row.names(raph.gene.counts.clean.subset))
raph.interaction.GO = topGO_wrapper(geneScores = GOscores.raph.interaction,
                                   geneScoresDE = F,
                                   geneScoresDirection = NA,
                                   GOmapping = GOmapping.raph,
                                   algorithm = "weight01",
                                   statistic = "fisher",
                                   nodeSize = 10,
                                   discretisedDE = F,
                                   p = 0.05)
raph.interaction.GO$consolidated_result
#26 GO terms: many epigenetic "histone binding" and "histone methylation", "chromatin binding", "regulation of histone modification"
#several response to sugar (fructose, glucose, sucrose), many developmental

#for the interaction terms, we also plot the output to understand what exactly is going on
#get the 12 terms with lowest adjusted p value
degs.raph.interaction.signif = subset(degs.raph.interaction, padj<0.05)
interestgenes = row.names(degs.raph.interaction.signif)[order(degs.raph.interaction.signif$log2FoldChange,decreasing = T)[1:12]]
#for each gene, extract the counts for plotting and label by gene
for(i in 1:length(interestgenes)){
  #if first element, instantiate frame, otherwise rbind to frame
  if(i==1) {
    raph.intplotdata = plotCounts(dds.gene.deg, gene=interestgenes[i], intgroup=c("domesticated","treatment"), returnData = T)
    raph.intplotdata$gene = interestgenes[i]
  } else {
    addrows = data.frame(plotCounts(dds.gene.deg, gene=interestgenes[i], intgroup=c("domesticated","treatment"), returnData = T), 
                         gene = interestgenes[i])
    raph.intplotdata = rbind(raph.intplotdata, addrows)
    }
}

#reorder for aesthetic reasons
raph.intplotdata$treatment = fct_rev(raph.intplotdata$treatment)
#plot with ggplot facet wrap
raph.intplot= ggplot(raph.intplotdata, aes(x = treatment, y = count)) +
  stat_summary(aes(group = domesticated), fun = mean, geom = "path") +
  stat_summary(aes(color = domesticated), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(color = domesticated), fun = mean, geom = "point", size = 4) +
  geom_point(aes(color = domesticated), size = 2) +
  scale_color_manual(values = brewer.pal(6,"Set3")[c(5,6)]) +
  facet_wrap(~gene, scales = "free")

#save plot
ggsave(raph.intplot, 
       filename = "raphanus_interaction_deg_plots.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  40, height = 25, units = "cm")


#######################

#now compare raphanistrum vs other wild
#subset to remove all but raphanistrum sativus and its wild ancestor
metadata.raph.wilds = subset(metadata.raph, species != "Raphanus sativus")
raph.gene.counts.clean.wilds = raph.gene.counts.clean[,as.character(metadata.raph.wilds$sample)]
#check that metadata and count matrices conform
table(metadata.raph.wilds$sample %in% colnames(raph.gene.counts.clean.wilds))
raph.gene.counts.clean.wilds = raph.gene.counts.clean.wilds[,as.character(metadata.raph.wilds$sample)]
#add column to check whether ancestor of domesticated or not
metadata.raph.wilds$wild.ancestor = (metadata.raph.wilds$species=="Raphanus raphanistrum")

#now run model with our variables of interest
dds.gene.wilds = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean.wilds,
                                        colData = metadata.raph.wilds,
                                        design = as.formula(~treatment+wild.ancestor+treatment*wild.ancestor))
dds.gene.deg.wilds = DESeq(dds.gene.wilds, fitType = "parametric", betaPrior = FALSE)

#results: interaction
degs.raph.wilds.interaction = results(dds.gene.deg.wilds, 
                                name="treatmentWheat.wild.ancestorTRUE",     
                                alpha = 0.05,
                                lfcThreshold = log2(1))
summary(degs.raph.wilds.interaction) #107 degs
degs.raph.interaction.ids = rownames(subset(degs.raph.interaction, padj<=0.05))
#5 of these degs are shared with parental degs
table(degs.raph.interaction.ids%in%parental.degs.ids)
# #also check for GO terms 
# raph.interaction.GO = topGO_wrapper(geneScores = GOscores.raph.interaction,
#                                     geneScoresDE = F,
#                                     geneScoresDirection = NA,
#                                     GOmapping = GOmapping.raph,
#                                     algorithm = "weight01",
#                                     statistic = "fisher",
#                                     nodeSize = 10,
#                                     discretisedDE = F,
#                                     p = 0.05)
# raph.interaction.GO$consolidated_result


#for the interaction terms, we also plot the output to understand what exactly is going on
#get the 12 terms with lowest adjusted p value
degs.raph.wilds.interaction.signif = subset(degs.raph.interaction, padj<0.05)
interestgeneswild = row.names(degs.raph.wilds.interaction.signif)[order(degs.raph.wilds.interaction.signif$log2FoldChange,decreasing = T)[1:12]]
#for each gene, extract the counts for plotting and label by gene
for(i in 1:length(interestgeneswild)){
  #if first element, instantiate frame, otherwise rbind to frame
  if(i==1) {
    raph.intplotdata = plotCounts(dds.gene.deg.wilds, gene=interestgeneswild[i], intgroup=c("wild.ancestor","treatment"), returnData = T)
    raph.intplotdata$gene = interestgeneswild[i]
  } else {
    addrows = data.frame(plotCounts(dds.gene.deg.wilds, gene=interestgeneswild[i], intgroup=c("wild.ancestor","treatment"), returnData = T), 
                         gene = interestgeneswild[i])
    raph.intplotdata = rbind(raph.intplotdata, addrows)
  }
}

#reorder for aesthetic reasons
raph.intplotdata$treatment = fct_rev(raph.intplotdata$treatment)
#plot with ggplot facet wrap
raph.intplot.wilds= ggplot(raph.intplotdata, aes(x = treatment, y = count)) +
  stat_summary(aes(group = wild.ancestor), fun = mean, geom = "path") +
  stat_summary(aes(color = wild.ancestor), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(color = wild.ancestor), fun = mean, geom = "point", size = 4) +
  geom_point(aes(color = wild.ancestor), size = 2) +
  scale_color_manual(values = brewer.pal(6,"Set3")[c(5,6)]) +
  facet_wrap(~gene, scales = "free") +
  labs(color = "Raphanistrum")

#save plot
ggsave(raph.intplot.wilds, 
       filename = "raphanus_interaction_deg_plots_wilds.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  40, height = 25, units = "cm")
