library(factoextra)
library(tidyverse)
library(cluster)
library(DESeq2)
library(ggpubr)

metadata.brass.wildrapa = subset(metadata.brass.clean,species=="Brassica rapa" & domesticated=="Wild")
brass.gene.counts.clean.wildrapa = brass.gene.counts.clean[,as.character(metadata.brass.wildrapa$sample)]
metadata.brass.domrapa = subset(metadata.brass.clean,species=="Brassica rapa" & domesticated!="Wild")
brass.gene.counts.clean.domrapa = brass.gene.counts.clean[,as.character(metadata.brass.domrapa$sample)]

#now run proper model
dds.gene = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.domrapa,
                                  colData = metadata.brass.domrapa,
                                  design = as.formula(~treatment))
dds.gene.deg = DESeq(dds.gene, fitType = "parametric", betaPrior = FALSE)
degs.brassdom.treatment = results(dds.gene.deg, 
                                  name="treatment_Control_vs_Wheat",     
                                  alpha = 0.05,
                                  lfcThreshold = log2(1))

#now run proper model
dds.gene = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.wildrapa,
                                  colData = metadata.brass.wildrapa,
                                  design = as.formula(~treatment))
dds.gene.deg = DESeq(dds.gene, fitType = "parametric", betaPrior = FALSE)
degs.brasswild.treatment = results(dds.gene.deg, 
                                   name="treatment_Control_vs_Wheat",     
                                   alpha = 0.05,
                                   lfcThreshold = log2(1))

#compile fold changes
brass.FC.frame = cbind(abs(degs.brasswild.treatment$log2FoldChange), 
                       abs(degs.brassdom.treatment$log2FoldChange)) %>% 
  'rownames<-'(row.names(degs.brasswild.treatment)) %>%
  'colnames<-'(c("Wild","Dom")) %>%
  na.omit() %>% data.frame()
#or use signed rather than absolute LFCs
brass.FC.frame = cbind(degs.brasswild.treatment$log2FoldChange, 
                       degs.brassdom.treatment$log2FoldChange) %>% 
  'rownames<-'(row.names(degs.brasswild.treatment)) %>%
  'colnames<-'(c("Wild","Dom")) %>%
  na.omit() %>% data.frame()

brass.FC.frame = brass.FC.frame[degs.brass.interaction.ids,]

set.seed(123)
# function to compute total within-cluster sum of square 
wss <- function(k) {kmeans(brass.FC.frame, k, nstart = 10 )$tot.withinss}
# Compute and plot wss for k = 1 to k = 10
k.values <- 1:10
# extract wss for 2-10 clusters
wss_values <- map_dbl(k.values, wss)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# gap_stat <- clusGap(brass.FC.frame, FUN = kmeans, nstart = 25,
#                     K.max = 10, B = 10)

# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(brass.FC.frame, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(brass.FC.frame))
  mean(ss[, 3])
}
# Compute and plot wss for k = 2 to k = 15
k.values <- 2:8
# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)
plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")

brass.FC.frame.clust = brass.FC.frame
#based on these, use 5 clusters for directional data or 3 clusters for non-directional data
foo = kmeans(brass.FC.frame,4,nstart = 25)
brass.FC.frame.clust$cluster = as.factor(foo$cluster)

clustplot = ggplot(brass.FC.frame.clust,aes(x=Wild,y=Dom,colour=cluster))+
  geom_point() +
  xlab("Log2 fold change (Wild rapa)") +
  ylab("Log2 fold change (Domesticated rapa)") #+
  #xlim(0,10) + ylim(0,10)
  # theme(aspect.ratio = 0.15,
  #       axis.title.y = element_blank(),
  #       axis.text.y = element_blank())

clustplot
#clustplot.inter = clustplot
#clustplot.signed = clustplot

#max(subset(foo, cluster==1)$means)
#min(subset(foo, cluster==1)$means)

ggarrange(clustplot.signed,clustplot.inter,nrow = 1,labels = c("A","B"))
