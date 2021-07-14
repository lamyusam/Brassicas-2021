extradata = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info.csv")
metadata.raph.extra = mutate(metadata.raph, extra = subset(extradata, substr(Species,1,3)=="Rap")$Wild..Domesticated)

#plot
ggplot(ggpcadata, aes(x = PC1, y = PC2, color = species, shape = domesticated, label = sample)) +
  geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
  geom_text(vjust = -1) +
  xlab(paste0("PC",1,": ",signif(pca.out$percent.var[1]*100, 3),"%")) +
  ylab(paste0("PC",2,": ",signif(pca.out$percent.var[2]*100, 3),"%")) +
  theme_bw() +
  scale_color_manual(name = "Treatment",
                     values = brewer.pal(7, "Paired")) +
  scale_shape_manual(name = "Treatment",
                     values = c(8,15:20)) +
  theme(panel.grid = element_line(color = "grey95"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(face = "bold", size =12))

#extract highest-variance genes
raph.gene.counts.clean.highvars = raph.gene.counts.clean[order(rowVars(as.matrix(raph.gene.counts.clean)),decreasing = T)[1:30000],] %>%
  select(-c("A110"))

# log-transform with a pseudocount
pca.counts = log2(raph.gene.counts.clean.highvars+1)
#create pca object
data.pca = prcomp(t(pca.counts), scale. = F)
#extract PC data
percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
pca.out = list(values = data.frame(data.pca$x),
               percent.var = percent.var)
#connect to phenotypic data
ggpcadata = pca.out$values %>%
  rownames_to_column(var = "sample") %>%
  left_join(metadata.raph.extra,
            by = "sample")

#plot
ggplot(ggpcadata, aes(x = PC1, y = PC2, color = species, shape = extra, label = sample)) +
  geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
  geom_text(vjust = -1) +
  xlab(paste0("PC",1,": ",signif(pca.out$percent.var[1]*100, 3),"%")) +
  ylab(paste0("PC",2,": ",signif(pca.out$percent.var[2]*100, 3),"%")) +
  theme_bw() +
  scale_color_manual(name = "Species",
                     values = brewer.pal(7, "Paired")) +
  labs(shape = "GbS category")
  theme(panel.grid = element_line(color = "grey95"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(face = "bold", size =12))

expression.heatmap(countdata = raph.gene.counts.clean.highvars,
                   data.phenotype = metadata.raph,
                   labels = c("species"),
                   pass_on = F,
                   ID_var = "sample")



####


extradata[grep("caudatus",extradata$Species),] %>%
  subset(!(RNAseq.sample.name %in% c("A15","A68","A36","A40")))

subset(extradata, Species == "Raphanus sativus")
  
