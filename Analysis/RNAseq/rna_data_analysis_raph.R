#### setup ####
# get libraries
basic_libraries <- c("DESeq2",
                     "tidyverse",
                     "WGCNA",
                     "RColorBrewer",
                     "forcats")

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
  dplyr::select(c("RNAseq.sample.name","Species","Parental.effects.status","Environment","Domesticated")) %>%
  'colnames<-'(c("sample","species","parental.effects","treatment","domesticated")) %>%
  mutate_all(as.factor) %>%
  mutate(treatment = fct_relevel(treatment, c("Wheat","Control")),
         domesticated = fct_relevel(domesticated, c("Wild","Cultivated")))
  

#subset by genus
metadata.brass = subset(data.meta, substr(species,1,3)=="Bra")
metadata.raph = subset(data.meta, substr(species,1,3)=="Rap")

#save for mark
write.csv(metadata.raph, file = "Data/RNAseq/metadata_raph.csv", row.names = F)

#pull rnaseq data
raph.gene.counts = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/raph.gene.counts.csv", row.names = 1)

#### filtering and qc ####
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
ggplot(ggpcadata, aes(x = PC1, y = PC2, color = species, shape = treatment, label = sample)) +
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

#### parental effects ####
#next step is to check for genes with parental effects
dds.parental = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean,
                                      colData = metadata.raph,
                                      design = as.formula(~parental.effects))
dds.parental.deg = DESeq(dds.parental, fitType = "parametric", betaPrior = FALSE)
parental.degs = results(dds.parental.deg,alpha = 0.1)
#very few genes have even marginal evidence of parental effects
print(paste0("Number of genes with parental effects at p<0.1: ",length(which(parental.degs$padj<0.1)),"/",nrow(parental.degs)))
parental.degs.ids = rownames(subset(parental.degs, padj<=0.05))
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


#### raphanistrum vs sativus ####
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
boxplot(log10(assays(dds.gene.deg)[["cooks"]]), range=0, las=2)
boxplot(log10(assays(dds.gene.deg)[["counts"]]), range=0, las=2)

#Also plot the dispersion extimates to make sure they look fine
plotDispEsts(dds.gene.deg)

#check effect of a parametric rather than local fit
parametric.gene.deg = DESeq(dds.gene, fitType = "parametric", betaPrior = FALSE)
plotDispEsts(parametric.gene.deg) 

# parametric seems to follow the data better so we'll use that going forward
dds.gene.deg = parametric.gene.deg

resultsNames(dds.gene.deg)

#results: wheat vs control
degs.raph.treatment = results(dds.gene.deg, 
                              name="treatment_Control_vs_Wheat",     
                              alpha = 0.05,
                              lfcThreshold = log2(1))
degs.raph.treatment.Nup = nrow(subset(degs.raph.treatment, padj<=0.05 & log2FoldChange>0)) #35 up in control
degs.raph.treatment.Ndown = nrow(subset(degs.raph.treatment, padj<=0.05 & log2FoldChange<0)) #47 up in wheat
degs.raph.treatment.ids = rownames(subset(degs.raph.treatment, padj<=0.05))
#2 of these degs shared with parental degs
table(degs.raph.treatment.ids%in%parental.degs.ids)
#also check for  GO terms 
raph.treatment.GO.up = topGO_wrapper(geneScores = degs.raph.treatment,
                         geneScoresDE = T,
                         geneScoresDirection = "Up",
                         GOmapping = GOmapping.raph,
                         algorithm = "weight01",
                         statistic = "fisher",
                         nodeSize = 10,
                         discretisedDE = T,
                         p = 0.05)
write.csv(raph.treatment.GO.up$consolidated_result, 
          file = "Analysis/RNAseq/Tables/raphanistrum_sativus_GO_controlbias.csv", row.names = FALSE)
#38 GO terms up
raph.treatment.GO.down = topGO_wrapper(geneScores = degs.raph.treatment,
                                     geneScoresDE = T,
                                     geneScoresDirection = "Down",
                                     GOmapping = GOmapping.raph,
                                     algorithm = "weight01",
                                     statistic = "fisher",
                                     nodeSize = 10,
                                     discretisedDE = T,
                                     p = 0.05)
write.csv(raph.treatment.GO.down$consolidated_result, 
          file = "Analysis/RNAseq/Tables/raphanistrum_sativus_GO_wheatbias.csv", row.names = FALSE)
#28 GO terms down

#results: cultivated vs wild
degs.raph.cultivated = results(dds.gene.deg, 
                               name="domesticated_Cultivated_vs_Wild",     
                               alpha = 0.05,
                               lfcThreshold = log2(1))
summary(degs.raph.cultivated) #~300 degs
degs.raph.cultivated.Nup = nrow(subset(degs.raph.cultivated, padj<=0.05 & log2FoldChange>0)) #75 up in cultivar
degs.raph.cultivated.Ndown = nrow(subset(degs.raph.cultivated, padj<=0.05 & log2FoldChange<0)) #203 up in wild
degs.raph.cultivated.ids = rownames(subset(degs.raph.cultivated, padj<=0.05))
#4 of these degs is shared with parental degs
table(degs.raph.cultivated.ids%in%parental.degs.ids)
#also check for deg GO terms 
raph.cultivated.GO.up = topGO_wrapper(geneScores = degs.raph.cultivated,
                                  geneScoresDE = T,
                                  geneScoresDirection = "Up",
                                  GOmapping = GOmapping.raph,
                                  algorithm = "weight01",
                                  statistic = "fisher",
                                  nodeSize = 10,
                                  discretisedDE = T,
                                  p = 0.05)
write.csv(raph.cultivated.GO.up$consolidated_result, 
          file = "Analysis/RNAseq/Tables/raphanistrum_sativus_GO_cultivatedbias.csv", row.names = FALSE)
#26 GO terms up
raph.cultivated.GO.down = topGO_wrapper(geneScores = degs.raph.cultivated,
                                      geneScoresDE = T,
                                      geneScoresDirection = "Down",
                                      GOmapping = GOmapping.raph,
                                      algorithm = "weight01",
                                      statistic = "fisher",
                                      nodeSize = 10,
                                      discretisedDE = T,
                                      p = 0.05)
write.csv(raph.cultivated.GO.down$consolidated_result, 
          file = "Analysis/RNAseq/Tables/raphanistrum_sativus_GO_wildbias.csv", row.names = FALSE)
#28 down

#results: interaction
degs.raph.interaction = results(dds.gene.deg, 
                                name="treatmentControl.domesticatedCultivated",     
                                alpha = 0.05,
                                lfcThreshold = log2(1))
summary(degs.raph.interaction) #173 degs
degs.raph.interaction.ids = rownames(subset(degs.raph.interaction, padj<=0.05))
degs.raph.interaction.N = length(degs.raph.interaction.ids)
#4 of these degs are shared with parental degs
table(degs.raph.interaction.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
raph.interaction.GO = topGO_wrapper(geneScores = degs.raph.interaction,
                                   geneScoresDE = T,
                                   geneScoresDirection = NA,
                                   GOmapping = GOmapping.raph,
                                   algorithm = "weight01",
                                   statistic = "fisher",
                                   nodeSize = 10,
                                   discretisedDE = T,
                                   p = 0.05)
raph.interaction.GO$consolidated_result
#57 GO terms

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
       filename = "raphanistrum_sativus_interaction_deg_plots.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  40, height = 25, units = "cm")

numdegs = c(degs.raph.cultivated.Nup, degs.raph.cultivated.Ndown, 
            degs.raph.treatment.Nup, degs.raph.treatment.Ndown, degs.raph.interaction.N)
numGO = c(nrow(raph.cultivated.GO.up$consolidated_result), nrow(raph.cultivated.GO.down$consolidated_result),
          nrow(raph.treatment.GO.up$consolidated_result), nrow(raph.treatment.GO.down$consolidated_result),
          nrow(raph.interaction.GO$consolidated_result))

raphanistrum.sativus.output = data.frame(DEGs=numdegs, GO_terms=numGO, 
                                          row.names = c("Domesticated_bias","Wild_bias","Unstressed_bias","Stressed_bias","Interaction"))

write.csv(raphanistrum.sativus.output, file = "Analysis/RNAseq/Tables/raphanistrum_sativus_summary.csv")

#### sativus + raphanistrum interaction norm analysis ####
#one thing we'd like to know: do genes generally gain or lose plasticity after domestication?
#get just sativus data
metadata.raph.sativus = subset(metadata.raph.subset, species == "Raphanus sativus")
raph.gene.counts.clean.sativus = raph.gene.counts.clean[,as.character(metadata.raph.sativus$sample)]
#get just raphanistrum data
metadata.raph.raphanistrum = subset(metadata.raph.subset, species == "Raphanus raphanistrum")
raph.gene.counts.clean.raphanistrum = raph.gene.counts.clean[,as.character(metadata.raph.raphanistrum$sample)]
#generate fold changes for sativus
dds.gene.sativus = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean.sativus,
                                             colData = metadata.raph.sativus,
                                             design = as.formula(~treatment))
dds.gene.deg.sativus = DESeq(dds.gene.sativus, fitType = "parametric", betaPrior = FALSE)
degs.sativus = results(dds.gene.deg.sativus,
                          name="treatment_Control_vs_Wheat",
                          alpha = 0.05,
                          lfcThreshold = log2(1))
#generate fold changes for raphanistrum
dds.gene.raphanistrum = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean.raphanistrum,
                                               colData = metadata.raph.raphanistrum,
                                               design = as.formula(~treatment))
dds.gene.deg.raphanistrum = DESeq(dds.gene.raphanistrum, fitType = "parametric", betaPrior = FALSE)
degs.raphanistrum = results(dds.gene.deg.raphanistrum,
                            name="treatment_Control_vs_Wheat",
                            alpha = 0.05,
                            lfcThreshold = log2(1))
#compile foldchange data for the significant interaction terms
foldchanges.subset = data.frame(raphanistrum = degs.raphanistrum[degs.raph.interaction.ids,"log2FoldChange"],
                               sativus = degs.sativus[degs.raph.interaction.ids,"log2FoldChange"], 
                               row.names = degs.raph.interaction.ids)
#in general, sativus genes are much more plastic
t.test(x = abs(foldchanges.subset$raphanistrum), y = abs(foldchanges.subset$sativus), paired = TRUE)
wilcox.test(x = abs(foldchanges.subset$raphanistrum), y = abs(foldchanges.subset$sativus), paired = TRUE)
mean(abs(foldchanges.subset$raphanistrum), na.rm = T)
mean(abs(foldchanges.subset$sativus), na.rm = T)
#plot
gg.foldchanges.subset = ggplot(data = reshape2::melt(foldchanges.subset), aes(y = abs(value), x =variable)) +
  geom_boxplot() + 
  #geom_point() +
  labs(x = "group", y = "absolute log2 fold change")
#save plot
ggsave(gg.foldchanges.subset, 
       filename = "raphanistrum_subset_foldchanges_plot.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  30, height = 25, units = "cm")

# another thing we'd like to know: have genes gained or lost plasticity in domestication, 
# and is the direction of that plasticity the same or opposite following domestication?

#pull the relevant data
foldchanges.subset = mutate(foldchanges.subset, 
                            direction=ifelse((sign(raphanistrum)==sign(sativus)),"Equal","Opposite"),
                            magnitude=ifelse((abs(raphanistrum)>abs(sativus)),"Decrease","Increase"))
table.foldchanges = table(select(foldchanges.subset,c("direction","magnitude")))
#distribution of opposite and equal changes is unrelated to magnitudes of changes:
chisq.test(table.foldchanges)
mosaicplot(table.foldchanges)

#more sophisticated
foldchanges.subset$wildDirection = ifelse(foldchanges.subset$wild>0.5, "up", 
                                          ifelse(foldchanges.subset$wild<(-0.5), "down","neutral"))
foldchanges.subset$cultivatedDirection = ifelse(foldchanges.subset$cultivated>0.5, "up", 
                                                ifelse(foldchanges.subset$cultivated<(-0.5), "down","neutral"))
table.foldchanges.comp = table(select(foldchanges.subset,c("wildDirection","cultivatedDirection")))
#test
chisq.test(table.foldchanges.comp)
mosaicplot(table.foldchanges.comp)

#oppposite direction is significantly more common than same direction
chisq.test(table(select(foldchanges.subset,c("direction"))))
#to confirm, increase in plasticity is significantly more common than decrease
chisq.test(table(select(foldchanges.subset,c("magnitude"))))

# just out of interest, it would be interesting to know whether these overall patterns also hold when looking at 
# the entire set of genes, not just the ones ID'd by DESeq2
#compile foldchange data for all terms
foldchanges.subset.all = data.frame(raphanistrum = degs.raphanistrum[,"log2FoldChange"],
                                    sativus = degs.sativus[,"log2FoldChange"], 
                                    row.names = row.names(degs.raphanistrum))
#yes, sativus genes are also more plastic overall
t.test(x = abs(foldchanges.subset.all$raphanistrum), y = abs(foldchanges.subset.all$sativus), paired = TRUE)
wilcox.test(x = abs(foldchanges.subset.all$raphanistrum), y = abs(foldchanges.subset.all$sativus), paired = TRUE)
mean(abs(foldchanges.subset.all$raphanistrum), na.rm = T)
mean(abs(foldchanges.subset.all$sativus), na.rm = T)
foldchanges.subset.all = mutate(foldchanges.subset.all, 
                            direction=ifelse((sign(raphanistrum)==sign(sativus)),"Equal","Opposite"),
                            magnitude=ifelse((abs(raphanistrum)>abs(sativus)),"Decrease","Increase"))
#intriguingly, when looking across all genes, same direction is more common than opposite direction
chisq.test(table(select(foldchanges.subset.all,c("direction"))))
#however, increase in plasticity is definitely more common even across all genes
chisq.test(table(select(foldchanges.subset.all,c("magnitude"))))

#### raphanistrum vs wilds ####

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
# Check for outliers: none are apparent
boxplot(log10(assays(dds.gene.deg.wilds)[["cooks"]]), range=0, las=2)
boxplot(log10(assays(dds.gene.deg.wilds)[["counts"]]), range=0, las=2)

#results: wheat vs control
degs.raph.wilds.treatment = results(dds.gene.deg.wilds, 
                              name="treatment_Control_vs_Wheat",     
                              alpha = 0.05,
                              lfcThreshold = log2(1))
degs.raph.wilds.treatment.Nup = nrow(subset(degs.raph.wilds.treatment, padj<=0.05 & log2FoldChange>0)) #37 up in control
degs.raph.wilds.treatment.Ndown = nrow(subset(degs.raph.wilds.treatment, padj<=0.05 & log2FoldChange<0)) #46 up in wheat
degs.raph.wilds.treatment.ids = rownames(subset(degs.raph.wilds.treatment, padj<=0.05))
#3 of these degs shared with parental degs
table(degs.raph.wilds.treatment.ids%in%parental.degs.ids)
#also check for  GO terms 
raph.wilds.treatment.GO.up = topGO_wrapper(geneScores = degs.raph.wilds.treatment,
                                     geneScoresDE = T,
                                     geneScoresDirection = "Up",
                                     GOmapping = GOmapping.raph,
                                     algorithm = "weight01",
                                     statistic = "fisher",
                                     nodeSize = 10,
                                     discretisedDE = T,
                                     p = 0.05)
write.csv(raph.wilds.treatment.GO.up$consolidated_result, 
          file = "Analysis/RNAseq/Tables/raphanus_wilds_GO_controlbias.csv", row.names = FALSE)
#44 GO terms up
raph.wilds.treatment.GO.down = topGO_wrapper(geneScores = degs.raph.wilds.treatment,
                                       geneScoresDE = T,
                                       geneScoresDirection = "Down",
                                       GOmapping = GOmapping.raph,
                                       algorithm = "weight01",
                                       statistic = "fisher",
                                       nodeSize = 10,
                                       discretisedDE = T,
                                       p = 0.05)
write.csv(raph.wilds.treatment.GO.down$consolidated_result, 
          file = "Analysis/RNAseq/Tables/raphanus_wilds_GO_wheatbias.csv", row.names = FALSE)
#35 GO terms down

#results: cultivated vs wild
degs.raph.wilds.cultivated = results(dds.gene.deg.wilds, 
                               name="wild.ancestorTRUE",     
                               alpha = 0.05,
                               lfcThreshold = log2(1))
summary(degs.raph.wilds.cultivated) #~2000 degs
degs.raph.wilds.cultivated.Nup = nrow(subset(degs.raph.wilds.cultivated, padj<=0.05 & log2FoldChange>0)) #734 up in cultivar
degs.raph.wilds.cultivated.Ndown = nrow(subset(degs.raph.wilds.cultivated, padj<=0.05 & log2FoldChange<0)) #1139 up in wild
degs.raph.wilds.cultivated.ids = rownames(subset(degs.raph.wilds.cultivated, padj<=0.05))
#12 of these degs is shared with parental degs
table(degs.raph.wilds.cultivated.ids%in%parental.degs.ids)
#also check for deg GO terms 
raph.wilds.cultivated.GO.up = topGO_wrapper(geneScores = degs.raph.wilds.cultivated,
                                      geneScoresDE = T,
                                      geneScoresDirection = "Up",
                                      GOmapping = GOmapping.raph,
                                      algorithm = "weight01",
                                      statistic = "fisher",
                                      nodeSize = 10,
                                      discretisedDE = T,
                                      p = 0.05)
write.csv(raph.wilds.cultivated.GO.up$consolidated_result, 
          file = "Analysis/RNAseq/Tables/raphanus_wilds_GO_cultivatedbias.csv", row.names = FALSE)
#68 GO terms up
raph.wilds.cultivated.GO.down = topGO_wrapper(geneScores = degs.raph.wilds.cultivated,
                                        geneScoresDE = T,
                                        geneScoresDirection = "Down",
                                        GOmapping = GOmapping.raph,
                                        algorithm = "weight01",
                                        statistic = "fisher",
                                        nodeSize = 10,
                                        discretisedDE = T,
                                        p = 0.05)
write.csv(raph.wilds.cultivated.GO.down$consolidated_result, 
          file = "Analysis/RNAseq/Tables/raphanus_wilds_GO_wildbias.csv", row.names = FALSE)
#69 down

#results: interaction
degs.raph.wilds.interaction = results(dds.gene.deg.wilds, 
                                name="treatmentControl.wild.ancestorTRUE",     
                                alpha = 0.05,
                                lfcThreshold = log2(1))
summary(degs.raph.wilds.interaction) #100 degs
degs.raph.wilds.interaction.ids = rownames(subset(degs.raph.wilds.interaction, padj<=0.05))
degs.raph.wilds.interaction.N = length(degs.raph.wilds.interaction.ids)
#2 of these degs are shared with parental degs
table(degs.raph.wilds.interaction.ids%in%parental.degs.ids)
#also check for GO terms
raph.wilds.interaction.GO = topGO_wrapper(geneScores = degs.raph.wilds.interaction,
                                    geneScoresDE = T,
                                    geneScoresDirection = NA,
                                    GOmapping = GOmapping.raph,
                                    algorithm = "weight01",
                                    statistic = "fisher",
                                    nodeSize = 10,
                                    discretisedDE = T,
                                    p = 0.05)
nrow(raph.wilds.interaction.GO$consolidated_result)

#for the interaction terms, we also plot the output to understand what exactly is going on
#get the 12 terms with lowest adjusted p value
degs.raph.wilds.interaction.signif = subset(degs.raph.wilds.interaction, padj<0.05)
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

#data for this need to be instantiated above
wildsnumdegs = c(degs.raph.wilds.cultivated.Nup, degs.raph.wilds.cultivated.Ndown,
            degs.raph.wilds.treatment.Nup, degs.raph.wilds.treatment.Ndown, degs.raph.wilds.interaction.N)
wildsnumGO = c(nrow(raph.wilds.cultivated.GO.up$consolidated_result), nrow(raph.wilds.cultivated.GO.down$consolidated_result),
          nrow(raph.wilds.treatment.GO.up$consolidated_result), nrow(raph.wilds.treatment.GO.down$consolidated_result),
          nrow(raph.wilds.interaction.GO$consolidated_result))

raph.wilds.output = data.frame(DEGs=wildsnumdegs, GO_terms=wildsnumGO,
                               row.names = c("Domesticated_bias","Wild_bias","Unstressed_bias","Stressed_bias","Interaction"))

write.csv(raph.wilds.output, file = "Analysis/RNAseq/Tables/raph_wilds_summary.csv")


#### wilds interaction norm analysis ####

#first compare treatments for raphanus and for wilds separately
#for each interaction gene
foo = row.names(degs.raph.wilds.interaction.signif)[1]
foo = "gene-LOC108805611"

#get fold changes in raphanistrum
#extract norm. counts for each combination
metadata.raph.raphanistrum = subset(metadata.raph.wilds, species == "Raphanus raphanistrum")
metadata.raph.otherwilds = subset(metadata.raph.wilds, species != "Raphanus raphanistrum")
# raphwheatcounts = counts(dds.gene.deg.wilds,normalized=T)[,as.character(subset(metadata.raph.raphanistrum, treatment == "Wheat")$sample)]
# raphcontrolcounts = counts(dds.gene.deg.wilds,normalized=T)[,as.character(subset(metadata.raph.raphanistrum, treatment == "Control")$sample)]
# otherwheatcounts = counts(dds.gene.deg.wilds,normalized=T)[,as.character(subset(metadata.raph.otherwilds, treatment == "Wheat")$sample)]
# othercontrolcounts = counts(dds.gene.deg.wilds,normalized=T)[,as.character(subset(metadata.raph.otherwilds, treatment == "Control")$sample)]
# #extract log2 fold changes for non-raphanistrum
# otherwilds.folds = log2(rowMeans(otherwheatcounts)/rowMeans(othercontrolcounts))
# raphanistrum.folds = log2(rowMeans(raphwheatcounts)/rowMeans(raphcontrolcounts))

raph.gene.counts.clean.otherwilds = raph.gene.counts.clean[,as.character(metadata.raph.otherwilds$sample)]
dds.gene.otherwilds = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean.otherwilds,
                                        colData = metadata.raph.otherwilds,
                                        design = as.formula(~treatment))
dds.gene.deg.otherwilds = DESeq(dds.gene.otherwilds, fitType = "parametric", betaPrior = FALSE)
degs.otherwilds = results(dds.gene.deg.otherwilds,
                          name="treatment_Control_vs_Wheat",
                          alpha = 0.05,
                          lfcThreshold = log2(1))

raph.gene.counts.clean.raphanistrum = raph.gene.counts.clean[,as.character(metadata.raph.raphanistrum$sample)]
dds.gene.raphanistrum = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean.raphanistrum,
                                             colData = metadata.raph.raphanistrum,
                                             design = as.formula(~treatment))
dds.gene.deg.raphanistrum = DESeq(dds.gene.raphanistrum, fitType = "parametric", betaPrior = FALSE)
degs.raphanistrum = results(dds.gene.deg.raphanistrum,
                          name="treatment_Control_vs_Wheat",
                          alpha = 0.05,
                          lfcThreshold = log2(1))

foldchanges.wilds = data.frame(raphanistrum = degs.raphanistrum[degs.raph.wilds.interaction.ids,"log2FoldChange"],
                               others = degs.otherwilds[degs.raph.wilds.interaction.ids,"log2FoldChange"], row.names = degs.raph.wilds.interaction.ids)

t.test(x = abs(foldchanges.wilds$raphanistrum), y = abs(foldchanges.wilds$others), paired = TRUE)
wilcox.test(x = abs(foldchanges.wilds$raphanistrum), y = abs(foldchanges.wilds$others), paired = TRUE)
mean(abs(foldchanges.wilds$raphanistrum), na.rm = T)
mean(abs(foldchanges.wilds$others), na.rm = T)
#plot
gg.foldchanges.wilds = ggplot(data = reshape2::melt(foldchanges.wilds), aes(y = abs(value), x =variable)) +
  geom_boxplot(range = 0) + 
  #geom_point() +
  labs(x = "group", y = "absolute log2 fold change")
gg.foldchanges.wilds
#save plot
ggsave(gg.foldchanges.wilds, 
       filename = "raphanistrum_wilds_foldchanges_plot.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  30, height = 25, units = "cm")

# fancier but uglier plot
# ggplot(reshape2::melt(as.matrix(foldchanges.wilds)), aes(y = abs(value))) +
#   geom_boxplot(aes(x = rep(c(-3, 3), each = nrow(foldchanges.wilds)), group = Var2), fill = 'steelblue') +
#   geom_point(aes(x = rep(c(-1, 1), each = nrow(foldchanges.wilds))), size = 5) +
#   geom_line(aes(x = rep(c(-1, 1), each = nrow(foldchanges.wilds)), group = Var1))

#also check whether raphanistrum genes generally have same directionality as other wilds
#pull the relevant data
table(select(foldchanges.wilds, c("raphanistrumDirection","otherWildsDirection")))
foldchanges.wilds = mutate(foldchanges.wilds, 
                            direction=ifelse((sign(raphanistrum)==sign(others)),"Equal","Opposite"),
                            magnitude=ifelse((abs(raphanistrum)<abs(others)),"Decrease","Increase"))
table.foldchanges = table(select(foldchanges.wilds,c("direction","magnitude")))
#distribution of opposite and equal changes is unrelated to magnitudes of changes:
chisq.test(table.foldchanges)
mosaicplot(table.foldchanges)

#more sophisticated directions
foldchanges.wilds$raphanistrumDirection = ifelse(foldchanges.wilds$raphanistrum>0.5, "up", 
                                                 ifelse(foldchanges.wilds$raphanistrum<(-0.5), "down","neutral"))
foldchanges.wilds$otherWildsDirection = ifelse(foldchanges.wilds$others>0.5, "up", 
                                               ifelse(foldchanges.wilds$others<(-0.5), "down","neutral"))
table.foldchanges.comp = select(foldchanges.wilds, c("raphanistrumDirection","otherWildsDirection"))
#test
chisq.test(table.foldchanges.comp)
mosaicplot(table.foldchanges.comp)


#oppposite direction is significantly more common than same direction
chisq.test(table(select(foldchanges.wilds,c("direction"))))
#to confirm, raphanistrum genes are more likely to exhibit increased plasticity relative to other wilds
chisq.test(table(select(foldchanges.wilds,c("magnitude"))))

