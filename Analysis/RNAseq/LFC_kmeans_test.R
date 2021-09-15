library(factoextra)
library(tidyverse)
library(cluster)

#extra accessions that need to be dropped:
phylogeny.drop = c("BCR-WP2",
                   "BCR-WP3",
                   "BIC-WP6",
                   "BMO-WP1",
                   "BMO-WP2")

foo = data.frame(means = abs(na.omit(degs.brass.cultivated.stress$log2FoldChange)))
foo = data.frame(means = abs(na.omit(degs.raph.cultivated.stress$log2FoldChange)))
foo = data.frame(means = abs(na.omit(degs.brassdom.treatment$log2FoldChange)))
foo = data.frame(means = abs(na.omit(degs.brasswild.treatment$log2FoldChange)))

foo = data.frame(means = abs(na.omit(degs.brassdom.treatment[degs.brass.interaction.ids,]$log2FoldChange)))


set.seed(123)
# function to compute total within-cluster sum of square 
wss <- function(k) {kmeans(foo, k, nstart = 10 )$tot.withinss}
# Compute and plot wss for k = 1 to k = 10
k.values <- 1:10
# extract wss for 2-10 clusters
wss_values <- map_dbl(k.values, wss)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

gap_stat <- clusGap(foo, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 10)

# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(foo$means, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(foo$means))
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


#foo = data.frame(means = rowMeans(brass.gene.counts.clean))
#fviz_nbclust(foo, kmeans, method = "wss")
bar = kmeans(foo$means,2,nstart = 25)


foo$cluster = as.factor(bar$cluster)

raphplot = ggplot(foo,aes(x=means,y=1,colour=cluster))+
  geom_point() +
  xlab("Log2 fold change") +
  scale_x_continuous(labels = scales::comma) +
  theme(aspect.ratio = 0.15,
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

raphplot

max(subset(foo, cluster==1)$means)
min(subset(foo, cluster==1)$means)

library(ggpubr)
ggarrange(brassplot,raphplot,ncol = 1,labels = c("Brassica","Raphanus"))
