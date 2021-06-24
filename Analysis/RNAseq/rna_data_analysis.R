library("DESeq2")
library("tidyverse")
library("WGCNA")
library("RColorBrewer")

setwd("/home/benjamin/Documents/Brassicas_repo")

#import functions
source("Functions/expression_heatmap.R")
source("Functions/topGO_wrapper.R")

#pull relevant metadata
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

#pull rnaseq data
raph.gene.counts = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/raph.gene.counts.csv", row.names = 1)

#filter by expression
raph.gene.counts.clean = raph.gene.counts[(which(rowMeans(raph.gene.counts)>=1)),]
paste0(nrow(raph.gene.counts)-nrow(raph.gene.counts.clean),"/",nrow(raph.gene.counts)," Raphanus genes filtered due to very low expression.")

#an additional filter- remove genes that have too many 0 counts
paste0("Removing an additional ",nrow(raph.gene.counts.clean[rowSums(raph.gene.counts.clean >= 5) < 3,])," genes with many low counts.")
raph.gene.counts.clean = raph.gene.counts.clean[rowSums(raph.gene.counts.clean >= 5) >= 3,]

#check that metadata and count matrices conform
table(metadata.raph$sample %in% colnames(raph.gene.counts.clean))
raph.gene.counts.clean = raph.gene.counts.clean[,as.character(metadata.raph$sample)]

#some basic QC: heatmaps and pca
# heatMap = expression.heatmap(countdata = sample_n(raph.gene.counts.clean,5000), 
#                              data.phenotype = metadata.raph,
#                              labels = c("parental.effects",
#                                         "treatment",
#                                         "domesticated"),
#                              pass_on = F,
#                              ID_var = "sample")

# # log-transform with a pseudocount
# pca.counts = log2(raph.gene.counts.clean+1)
# #create pca object
# data.pca = prcomp(t(pca.counts))
# #extract PC data
# percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
# pca.out = list(values = data.frame(data.pca$x),
#                percent.var = percent.var)
# #connect to phenotypic data
# ggpcadata = pca.out$values %>% 
#   rownames_to_column(var = "sample") %>%
#   left_join(metadata.raph,
#             by = "sample")
# #plot
# ggplot(ggpcadata, aes(x = PC1, y = PC2, color = domesticated, shape = treatment, label = sample)) +
#   geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
#   #geom_text(vjust = -1) +
#   xlab(paste0("PC",1,": ",signif(pca.out$percent.var[1]*100, 3),"%")) +
#   ylab(paste0("PC",2,": ",signif(pca.out$percent.var[2]*100, 3),"%")) +
#   theme_bw() +
#   # scale_color_manual(name = "Treatment",
#   #                    values = brewer.pal(7, "Paired")) +
#   scale_shape_manual(name = "Treatment",
#                      values = c(8,15:20)) +
#   theme(panel.grid = element_line(color = "grey95"),
#         legend.title = element_text(face = "bold"),
#         axis.text.x = element_text(size = 11),
#         axis.text.y = element_text(size = 11),
#         axis.title = element_text(face = "bold", size =12))


#next step is to check for genes with parental effects and exclude these
dds.parental = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean,
                                      colData = metadata.raph,
                                      design = as.formula(~parental.effects))
dds.parental.deg = DESeq(dds.parental, fitType = "local", betaPrior = FALSE)
parental.degs = results(dds.parental.deg)
#very few genes have even a hint of parental effect, so let's just drop these
print(paste0("Number of genes with parental effects at p<0.1: ",length(which(parental.degs$padj<0.1)),"/",nrow(parental.degs)))
parental.degs.ids = rownames(subset(parental.degs, padj<=0.1))
#raph.gene.counts.clean = raph.gene.counts.clean[rownames(subset(parental.degs, padj>=0.1)),]
#also check for parental deg GO terms 
GOscores.parental = as.numeric(row.names(raph.gene.counts.clean)%in%parental.degs.ids) %>% 'names<-'(row.names(raph.gene.counts.clean))
parentGO = topGO_wrapper(geneScores = GOscores.parental,
                         geneScoresDE = F,
                         geneScoresDirection = NA,
                         GOmapping = GOmapping.raph,
                         algorithm = "weight01",
                         statistic = "fisher",
                         nodeSize = 10,
                         discretisedDE = F,
                         p = 0.05)
parentGO$consolidated_result
#only 12 GO terms, with no clear interesting trend

#now run proper model
dds.gene = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean,
                                  colData = metadata.raph,
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
# parametric seems to follow the data better so let's go with that
dds.gene.deg = parametric.gene.deg

resultsNames(dds.gene.deg)



#wheat vs control
degs.raph.treatment = results(dds.gene.deg, 
                              name="treatment_Wheat_vs_Control",     
                              alpha = 0.05,
                              lfcThreshold = log2(1))
summary(degs.raph.treatment) #126 DEGs
degs.raph.treatment.ids = rownames(subset(degs.raph.treatment, padj<=0.05))
#1 of these degs is shared with parental degs
table(degs.raph.treatment.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
#GOscores.raph.treatment = as.numeric(row.names(raph.gene.counts.clean)%in%degs.raph.treatment.ids) %>% 'names<-'(row.names(raph.gene.counts.clean))
#terms up:
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

#cultivated vs wild
degs.raph.cultivated = results(dds.gene.deg, 
                               name="domesticated_Cultivated_vs_Wild",     
                               alpha = 0.05,
                               lfcThreshold = log2(1))
summary(degs.raph.cultivated) #~1500 degs
degs.raph.cultivated.ids = rownames(subset(degs.raph.cultivated, padj<=0.05))
#12 of these degs is shared with parental degs
table(degs.raph.cultivated.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
#GOscores.raph.cultivated = as.numeric(row.names(raph.gene.counts.clean)%in%degs.raph.cultivated.ids) %>% 'names<-'(row.names(raph.gene.counts.clean))
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


#interaction
degs.raph.interaction = results(dds.gene.deg, 
                                name="treatmentWheat.domesticatedCultivated",     
                                alpha = 0.05,
                                lfcThreshold = log2(1))
summary(degs.raph.interaction) #96 degs
degs.raph.interaction.ids = rownames(subset(degs.raph.interaction, padj<=0.05))
#1 of these degs is shared with parental degs
table(degs.raph.interaction.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
GOscores.raph.interaction = as.numeric(row.names(raph.gene.counts.clean)%in%degs.raph.interaction.ids) %>% 'names<-'(row.names(raph.gene.counts.clean))
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
#get the 12 terms with lowest padj
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

ggsave(raph.intplot, 
       filename = "raphanus_interaction_deg_plots.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  40, height = 25, units = "cm")


####### wgcna

# Perform filtering- we apply a more stringent filtering process here, because WGCNA benefits from cleaner data
# Here we remove any samples that don't have at least 5 counts in at least ~90% of samples
#SHOULD REALLY DO THIS ON COMPLETE DATA BEFORE SUBSETTING, ASSUMING WE WISH TO RUN THROUGH BRASS AND RAPH IN THE SAME ANALYSIS 
raph.gene.count.clean.wgcna = raph.gene.counts.clean[rowSums(raph.gene.counts.clean > 5) > (ncol(raph.gene.counts.clean)*0.9),]

# We also normalize the data prior to subsetting
raph.gene.count.clean.wgcna = as.matrix(raph.gene.count.clean.wgcna)  %>% 
  varianceStabilizingTransformation() %>% limma::normalizeBetweenArrays()

#subset data
wildData = raph.gene.count.clean.wgcna[,colnames(raph.gene.count.clean.wgcna)%in%subset(metadata.raph,domesticated=="Wild")$sample]
domesticatedData = raph.gene.count.clean.wgcna[,colnames(raph.gene.count.clean.wgcna)%in%subset(metadata.raph,domesticated=="Cultivated")$sample]
nSets = 2;

# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Wild","Domesticated")
shortLabels = c("wild","domestic")
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(wildData)))
names(multiExpr[[1]]$data) = rownames(wildData)
rownames(multiExpr[[1]]$data) = colnames(wildData)

multiExpr[[2]] = list(data = as.data.frame(t(domesticatedData)))
names(multiExpr[[2]]$data) = rownames(domesticatedData)
rownames(multiExpr[[2]]$data) = colnames(domesticatedData)

exprSize = checkSets(multiExpr)

gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

#cut genes without enough reads 
if (!gsg$allOK){
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

#cluster samples by euclidean distance (separately in each set)
sampleTrees = list()
for (set in 1:nSets){
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets){
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)
}


# Choose the cut height for the data set
baseHeight = 180
# Adjust the cut height for the data set for the number of samples?
#cutHeights = c(baseHeight, baseHeight*exprSize$nSamples[2]/exprSize$nSamples[1], baseHeight*exprSize$nSamples[3]/exprSize$nSamples[1]);
cutHeights = c(baseHeight,baseHeight)
# Choose the "base" cut height for the female data set
# Re-plot the dendrograms including the cut lines
#pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
for(set in 1:nSets){
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h=cutHeights[set], col = "red");
}

#cut outlier samples
for (set in 1:nSets){
  # Find clusters cut by the line
  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set], minSize = 10)
  # Keep the largest one (labeled by the number 1)
  keep = (labels==1)
  print(labels)
  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}
# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize

#attach phenotypic data
allTraits = dplyr::select(metadata.raph, c("sample",
                                           "treatment",
                                           "domesticated"))

# Form a multi-set structure that will hold the traits.
Traits = vector(mode="list", length = nSets)
for (set in 1:nSets){
  setSamples = rownames(multiExpr[[set]]$data)
  traitRows = match(setSamples, allTraits$sample)
  Traits[[set]] = list(data = allTraits[traitRows,])
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1]
}
collectGarbage(); #clean up memory
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

# # Choose a set of soft-thresholding powers
# #higher powers reduce heterogeneity more 
# powers = c(seq(4,10,by=1), seq(12,20, by=2));
# # Initialize a list to hold the results of scale-free analysis
# powerTables = vector(mode = "list", length = nSets);
# # Call the network topology analysis function for each set in turn
# for (set in 1:nSets){
#   powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, networkType="signed",
#                                                      verbose = 2)[[2]])
# }


# # Plot the results:
# colors = c("black", "red")
# # Will plot these columns of the returned scale free analysis tables
# plotCols = c(2,5,6,7)
# colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
#              "Max connectivity");
# # Get the minima and maxima of the plotted points
# ylim = matrix(NA, nrow = 2, ncol = 4);
# for (set in 1:nSets){
#   for (col in 1:length(plotCols)){
#     ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
#     ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
#   }
# }
# # Plot the quantities in the chosen columns vs. the soft thresholding power
# sizeGrWindow(8, 6)
# #pdf(file = "/home/benjamin/Dropbox/Ben PhD/Chapter_2_Manuscript/figures/softpowerplots.pdf", width = 18/2.54, height = 14/2.54)
# par(mfcol = c(2,2));
# par(mar = c(2.2, 2.2, 2.2, 2.2))
# cex1 = 1;
# for (col in 1:length(plotCols)) for (set in 1:nSets)
# {
#   if (set==1){
#     plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
#          xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
#          main = colNames[col]);
#     addGrid();
#   }
#   if (col==1)
#   {
#     text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
#          labels=powers,cex=cex1,col=colors[set]);
#   } else
#     text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
#          labels=powers,cex=cex1,col=colors[set]);
#   if (col==1)
#   {
#     legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
#   } else
#     legend("topright", legend = setLabels, col = colors, pch = 20) ;
# }
# dev.off()


# beepr::beep(3)
# 
# #lookslike 12 or 14 will be okay- this doesn't get us to a scale-free fit of 0.9, but it's the point where the fit levels off
# softPower = 14; 
# # Initialize an appropriate array to hold the adjacencies
# adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# # Calculate adjacencies in each individual data set
# for (set in 1:nSets){
#   #adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower
#   adjacencies[set, , ] = adjacency(multiExpr[[set]]$data, power = softPower, type = "signed")
# }
# 
# # Scaling of Topological Overlap Matrices to make them comparable across sets
# 
# # Initialize an appropriate array to hold the TOMs
# TOM = array(0, dim = c(nSets, nGenes, nGenes));
# # Calculate TOMs in each individual data set
# for (set in 1:nSets){
#   TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])
# }

#calculating TOMs is very slow, so here we just load the output of a pre-performed analysis directly
load("Analysis/RNAseq/raph_TOM_subsample8000.Rdata")

#get TOM dissimilarities
dissTOM = 1-pmin(TOM[1,,])

##### Clustering and module identification
# Clustering
consTree = hclust(as.dist(dissTOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = dissTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)

# sizeGrWindow(8,6)
# plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)

###### merge highly co-expressed modules
# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");

plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.4, col = "red")

merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.4, verbose = 3)

# Numeric module labels
moduleLabels = merge$colors;
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs

sizeGrWindow(8, 6)
png(file = "Analysis/RNAseq/Images/raph_WGCNA_subsample_mergedDendrosAndColors.png", wi = 18/2.54, he = 12/2.54, units = "in", res = 800)
par(mfcol = c(1,1));
par(mar = c(3.2, 3.2 , 3.2, 3.2))
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

######

exprSize = checkSets(multiExpr);
nSets = exprSize$nSets;


# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
# # Calculate the correlations
# for(set in 1:nSets){
#   #convert to numeric first to allow correlations
#   Traits[[set]]$data = dplyr::select(mutate_at(Traits[[set]]$data, .vars = c("treatment"), .funs = function(x) as.numeric(as.factor(x))),"treatment")
#   moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
#   moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
# }

set = 2

pvalues=c()
for(i in 1:length(consMEs[[set]]$data)) {
  
  glim = glm((as.numeric(Traits[[set]]$data$treatment)-1) ~ consMEs[[set]]$data[,i], family = "binomial")
  pval = summary(glim)$coefficients[2,4]
  pvalues = c(pvalues,pval)
}


# Convert numerical labels to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");

# Plot the module-trait relationship table for set number 1
set=1
options(scipen = 1)
padjmatrix_1 = matrix(p.adjust(signif(moduleTraitPvalue[[set]], 1),method = "BH"), 
                      nrow = nrow(moduleTraitPvalue[[set]]), 
                      ncol = ncol(moduleTraitPvalue[[set]]))
textMatrix_1 = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                     signif(c(padjmatrix_1),2), ")", sep = "");
dim(textMatrix_1) = dim(moduleTraitCor[[set]])

ColorNames = stringr::str_sub(MEColorNames,3,)

padjmelt_1 = 
  padjmatrix_1 %>% 
  data.frame() %>% 
  'colnames<-'(names(Traits[[set]]$data)) %>%
  'rownames<-'(ColorNames) %>%
  rownames_to_column(var = "modcolor") %>%
  reshape2::melt(id.vars = "modcolor") 

padjmelt_1$modcolor = factor(padjmelt_1$modcolor,levels = ColorNames)

NiceColorNames = paste0(toupper(str_sub(padjmelt_1$modcolor,1,1)),str_sub(padjmelt_1$modcolor,2,))
NiceColorNames = paste0(1:dim(textMatrix_1)[1],": ", NiceColorNames)

SimpleNames = paste0("Module ",c(29:1))

GGheat_1 = ggplot(data = padjmelt_1, aes(x = variable, y = modcolor)) +
  geom_tile(aes(fill = value),color = "gray", size=.75, width=0.5, height = 1) +
  geom_text(aes(label=c(textMatrix_1)), 
            lineheight = 0.75, size = 3.5) +
  # geom_point(data = padjmelt_1[1:22,],aes(x = 0.45, y = 1:22),
  #            fill = c(as.character(padjmelt_1$modcolor[1:22])),
  #            color = "black",
  #            size = 11,
  #            shape = 22) +
  scale_x_discrete("Trait",
                   expand = c(0,0),
                   labels = c("Wheat vs Control"),
                   position = "top") +
  scale_y_discrete("Module",
                   expand = c(0,0),
                   labels = SimpleNames,
                   limits = rev(levels(padjmelt_1$modcolor))) +
  scale_fill_gradientn(colours = colorRampPalette(rev(c("#FFFFFF",brewer.pal(n = 9, name = "Reds")[1:5])),bias=6)(20),
                       breaks = c(0.0,0.05,1),
                       expand = c(0,0),
                       limits = c(0,1),
                       guide = guide_colourbar(barheight = 25,
                                               title = "p-value\n(adjusted)",
                                               title.vjust = 2,
                                               frame.colour = "black", 
                                               frame.linewidth = 1.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size =11, color = "grey40"),
        axis.title = element_text(face = "bold", size =12),
        axis.title.x = element_blank())

GGheat_1

ggsave(GGheat_1, filename = "Raph_subset_moduletraitcor_wild.png", device = "png", path ="Analysis/RNAseq/Images/",
       width = 10, height = 20, units = "cm")

# Plot the module-trait relationship table for set number 2
set=2
options(scipen = 1)
padjmatrix_2 = matrix(p.adjust(signif(moduleTraitPvalue[[set]], 1),method = "BH"), 
                      nrow = nrow(moduleTraitPvalue[[set]]), 
                      ncol = ncol(moduleTraitPvalue[[set]]))
textMatrix_2 = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                     signif(c(padjmatrix_2),2), ")", sep = "");
dim(textMatrix_2) = dim(moduleTraitCor[[set]])

ColorNames = stringr::str_sub(MEColorNames,3,)

padjmelt_2 = 
  padjmatrix_2 %>% 
  data.frame() %>% 
  'colnames<-'(names(Traits[[set]]$data)) %>%
  'rownames<-'(ColorNames) %>%
  rownames_to_column(var = "modcolor") %>%
  reshape2::melt(id.vars = "modcolor") 

padjmelt_2$modcolor = factor(padjmelt_2$modcolor,levels = ColorNames)

NiceColorNames = paste0(toupper(str_sub(padjmelt_2$modcolor,1,1)),str_sub(padjmelt_2$modcolor,2,))
NiceColorNames = paste0(1:dim(textMatrix_2)[1],": ", NiceColorNames)

GGheat_2 = ggplot(data = padjmelt_2, aes(x = variable, y = modcolor)) +
  geom_tile(aes(fill = value),color = "gray", size=.75, width=1, height = 1) +
  geom_text(aes(label=c(textMatrix_2)), 
            lineheight = 0.75, size = 3.5) +
  # geom_point(data = padjmelt_2[1:22,],aes(x = 0.45, y = 1:22),
  #            fill = c(as.character(padjmelt_2$modcolor[1:22])),
  #            color = "black",
  #            size = 11,
  #            shape = 22) +
  scale_x_discrete("Trait",
                   expand = c(0,0),
                   labels = c("Wheat vs Control"),
                   position = "top") +
  scale_y_discrete("Module",
                   expand = c(0,0),
                   labels = SimpleNames,
                   limits = rev(levels(padjmelt_1$modcolor))) +
  scale_fill_gradientn(colours = colorRampPalette(rev(c("#FFFFFF",brewer.pal(n = 9, name = "Reds")[1:5])),bias=6)(20),
                       breaks = c(0.0,0.05,1),
                       expand = c(0,0),
                       limits = c(0,1),
                       guide = guide_colourbar(barheight = 25,
                                               title = "p-value\n(adjusted)",
                                               title.vjust = 2,
                                               frame.colour = "black", 
                                               frame.linewidth = 1.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size =11, color = "grey40"),
        axis.title = element_text(face = "bold", size =12),
        axis.title.x = element_blank())

GGheat_2
  
ggsave(GGheat_2, filename = "Raph_subset_moduletraitcor_domesticated.png", device = "png", path ="Analysis/RNAseq/Images/",
       width = 10, height = 20, units = "cm")
