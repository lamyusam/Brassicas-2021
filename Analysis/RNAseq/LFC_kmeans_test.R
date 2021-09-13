library(factoextra)

#extra accessions that need to be dropped:
phylogeny.drop = c("BCR-WP2",
                   "BCR-WP3",
                   "BIC-WP6",
                   "BMO-WP1",
                   "BMO-WP2")

foo = data.frame(means = abs(na.omit(degs.brass.treatment$log2FoldChange)))
foo = data.frame(means = abs(na.omit(degs.raph.treatment$log2FoldChange)))


#foo = data.frame(means = rowMeans(brass.gene.counts.clean))
bar = kmeans(foo$means,5,nstart = 25)
fviz_nbclust(bar, kmeans, method = "wss")


foo$cluster = as.factor(bar$cluster)

raphplot = ggplot(foo,aes(x=means,y=1,colour=cluster))+
  geom_point() +
  xlab("Log2 fold change") +
  scale_x_continuous(labels = scales::comma) +
  theme(aspect.ratio = 0.15,
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

raphplot

max(subset(foo, cluster==5)$means)
min(subset(foo, cluster==1)$means)

library(ggpubr)
ggarrange(brassplot,raphplot,ncol = 1,labels = c("Brassica","Raphanus"))
