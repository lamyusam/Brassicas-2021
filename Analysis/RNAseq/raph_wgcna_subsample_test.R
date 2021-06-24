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


# Perform filtering- we apply a more stringent filtering process here, because WGCNA benefits from cleaner data
# Here we remove any samples that don't have at least 5 counts in at least ~90% of samples
#SHOULD REALLY DO THIS ON COMPLETE DATA BEFORE SUBSETTING, ASSUMING WE WISH TO RUN THROUGH BRASS AND RAPH IN THE SAME ANALYSIS 
#raph.gene.count.clean.wgcna = raph.gene.counts.clean[rowSums(raph.gene.counts.clean > 5) > (ncol(raph.gene.counts.clean)*0.9),]
raph.gene.count.clean.wgcna = raph.gene.counts.clean[rowSums(raph.gene.counts.clean > 10) > (ncol(raph.gene.counts.clean)*0.95),]
#for now, subsample for speed
raph.gene.count.clean.wgcna = slice_sample(raph.gene.count.clean.wgcna, n = 8000)

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
baseHeight = 110
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

# Choose a set of soft-thresholding powers
#higher powers reduce heterogeneity more
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets){
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, networkType="signed",
                                                     verbose = 2)[[2]])
}


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


beepr::beep(3)

#lookslike 12 or 14 will be okay- this doesn't get us to a scale-free fit of 0.9, but it's the point where the fit levels off
softPower = 14;
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets){
  #adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower
  adjacencies[set, , ] = adjacency(multiExpr[[set]]$data, power = softPower, type = "signed")
}

#save(adjacencies,file = "Analysis/RNAseq/raph_adjacencies.Rdata")

# load("Analysis/RNAseq/raph_adjacencies.Rdata")
# load("Analysis/RNAseq/raph_TOM.Rdata")
# beepr::beep(3)

# Scaling of Topological Overlap Matrices to make them comparable across sets

# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets){
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])
}

save(TOM,file = "Analysis/RNAseq/raph_TOM_subsample8000.Rdata")
