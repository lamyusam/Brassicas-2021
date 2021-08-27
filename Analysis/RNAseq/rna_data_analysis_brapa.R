#### setup ####

library("DESeq2")
library("tidyverse")
library("WGCNA")

setwd("/home/benjamin/Documents/Brassicas_repo")

#import functions
source("Functions/expression_heatmap.R")
source("Functions/topGO_wrapper.R")

#import GO mapping
load("Data/GO/brass_GOmapping.Rdata")

#pull relevant metadata
data.meta = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info.csv") %>%
  mutate(Domesticated = ifelse(Wild..Domesticated=="Wild","Wild","Cultivated")) %>%
  mutate(Environment = ifelse(Environment=="wheat competition","Wheat","Control"),
         Parental.effects.status = ifelse(Parental.effects.status == "'\"standardised\"","Standardised","Unstandardised")) %>%
  dplyr::select(c("RNAseq.sample.name","Species","Label","Parental.effects.status","Environment","Domesticated")) %>%
  #dplyr::select(c("RNAseq.sample.name","Species","Parental.effects.status","Environment","Wild..Domesticated")) %>%
  'colnames<-'(c("sample","species","label","parental.effects","treatment","domesticated")) %>%
  mutate_all(as.factor)

data.meta = mutate(data.meta, 
                   treatment = fct_relevel(treatment, c("Wheat","Control")),
                   domesticated = fct_relevel(domesticated, c("Wild","Cultivated")))


#subset by genus
metadata.brass = subset(data.meta, substr(species,1,3)=="Bra")
metadata.raph = subset(data.meta, substr(species,1,3)=="Rap")

#pull rnaseq data
brass.gene.counts = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/brapa.gene.counts.csv", row.names = 1)

#### filtering and qc ####
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
#this is unusual- we see a clear separation of samples along the first PC, but the source isn't clear
#replot but with species labels
ggplot(ggpcadata, aes(x = PC1, y = PC2, shape = domesticated, color = species, label = sample)) +
  geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
  geom_text(vjust = -1) +
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
#from this we can see that the big clstering is rapa vs others. For simplicity, let's consider dropping the 'weird' samples,
#and also their 'partners' where relevant
checkframe = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info.csv")
checksample = substr(subset(checkframe, RNAseq.sample.name == "A31")$Label,1,11)
brassica.outliers = subset(checkframe, substr(Label,1,11) == checksample)$RNAseq.sample.name %>%
  c("A112","A108","A86","A129")
metadata.brass.clean= subset(metadata.brass, !(sample %in% brassica.outliers))
brass.gene.counts.clean = brass.gene.counts.clean[,as.character(metadata.brass.clean$sample)]

#re-run PCA now that outliers have been excluded
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
  left_join(metadata.brass.clean,
            by = "sample")
#plot with parental effects
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
#replot but with species labels
ggplot(ggpcadata, aes(x = PC1, y = PC2, shape = domesticated, color = species, label = sample)) +
  geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
  geom_text(vjust = -1) +
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
#okay so now we have a new problem- the brapa samples clearly divide into two very distinct groups along the second PC
#all the 'outlier' samples are cultivated and most are standarised, but the main group is composed of a mix, so what's the source?
outgroup = subset(ggpcadata, PC2<(-150))$sample
subset(checkframe, RNAseq.sample.name %in% outgroup)
#okay, the 'outliers' comprise all and only the tricoloris samples. Mark says this is fine!

#### parental effects ####
#next step is to check for genes with parental effects and exclude these
#since parental effects are only available for brapa, we have to make sure we don't include other species
dds.parental = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean[,as.character(subset(metadata.brass.clean, 
                                                                                               (species == "Brassica rapa"))$sample)],
                                      colData = subset(metadata.brass.clean, species == "Brassica rapa"),
                                      design = as.formula(~parental.effects))
dds.parental.deg = DESeq(dds.parental, fitType = "local", betaPrior = FALSE)
parental.degs = results(dds.parental.deg)
#an unremarkable proportion of genes have parental effects
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
#38 GO terms, many biosynthetic

#### wild brapa vs cultivated brapa ####
#subset to remove all but raphanistrum sativus and its wild ancestor
metadata.brass.subset = subset(metadata.brass.clean, species == "Brassica rapa")
brass.gene.counts.clean.subset = brass.gene.counts.clean[,as.character(metadata.brass.subset$sample)]
#check that metadata and count matrices conform
table(metadata.brass.subset$sample %in% colnames(brass.gene.counts.clean.subset))
brass.gene.counts.clean.subset = brass.gene.counts.clean.subset[,as.character(metadata.brass.subset$sample)]

#now run proper model
dds.gene = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.subset,
                                  colData = metadata.brass.subset,
                                  design = as.formula(~treatment+domesticated+treatment*domesticated))

# Run the default analysis for DESeq2 and generate results table. 
#NA p-values are generated by 0 counts and outliers calculated by Cook's distance.
dds.gene.deg = DESeq(dds.gene, fitType = "local", betaPrior = FALSE)
# Check for outliers: none are apparent
print("Check for gene expression outliers")
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds.gene.deg)[["cooks"]]), range=0, las=2)
boxplot(log10(assays(dds.gene.deg)[["counts"]]), range=0, las=2)

#Also plot the dispersion extimates to make sure they look fine
plotDispEsts(dds.gene.deg)

#check effectof a parametric rather than local fit
parametric.gene.deg = DESeq(dds.gene, fitType = "parametric", betaPrior = FALSE)
plotDispEsts(parametric.gene.deg) 
# parametric seems to follow the data a little better so let's go with that
dds.gene.deg = parametric.gene.deg

#save for later
save(dds.gene.deg, file = "Analysis/RNAseq/Tables/brassica_rapa_deseq.R")

resultsNames(dds.gene.deg)

#results: wheat vs control
degs.brass.treatment = results(dds.gene.deg, 
                              name="treatment_Control_vs_Wheat",     
                              alpha = 0.05,
                              lfcThreshold = log2(1))
degs.brass.treatment.Nup = nrow(subset(degs.brass.treatment, padj<=0.05 & log2FoldChange>0)) #16 up in control
degs.brass.treatment.Ndown = nrow(subset(degs.brass.treatment, padj<=0.05 & log2FoldChange<0)) #22 up in wheat
degs.brass.treatment.ids = rownames(subset(degs.brass.treatment, padj<=0.05))
#1 of these degs shared with parental degs
table(degs.brass.treatment.ids%in%parental.degs.ids)
#also check for  GO terms 
brass.treatment.GO.up = topGO_wrapper(geneScores = degs.brass.treatment,
                                     geneScoresDE = T,
                                     geneScoresDirection = "Up",
                                     GOmapping = GOmapping.brass,
                                     algorithm = "weight01",
                                     statistic = "fisher",
                                     nodeSize = 10,
                                     discretisedDE = T,
                                     p = 0.05)
write.csv(brass.treatment.GO.up$consolidated_result, 
          file = "Analysis/RNAseq/Tables/brapas_GO_controlbias.csv", row.names = FALSE)
#14 GO terms up
brass.treatment.GO.down = topGO_wrapper(geneScores = degs.brass.treatment,
                                       geneScoresDE = T,
                                       geneScoresDirection = "Down",
                                       GOmapping = GOmapping.brass,
                                       algorithm = "weight01",
                                       statistic = "fisher",
                                       nodeSize = 10,
                                       discretisedDE = T,
                                       p = 0.05)
write.csv(brass.treatment.GO.down$consolidated_result, 
          file = "Analysis/RNAseq/Tables/brapas_GO_wheatbias.csv", row.names = FALSE)
#20 GO terms down

#results: cultivated vs wild
degs.brass.cultivated = results(dds.gene.deg, 
                               name="domesticated_Cultivated_vs_Wild",     
                               alpha = 0.05,
                               lfcThreshold = log2(1))
summary(degs.brass.cultivated) #~2000 degs
degs.brass.cultivated.Nup = nrow(subset(degs.brass.cultivated, padj<=0.05 & log2FoldChange>0)) #1393 up in cultivar
degs.brass.cultivated.Ndown = nrow(subset(degs.brass.cultivated, padj<=0.05 & log2FoldChange<0)) #771 up in wild
degs.brass.cultivated.ids = rownames(subset(degs.brass.cultivated, padj<=0.05))
#41 of these degs are shared with parental degs
table(degs.brass.cultivated.ids%in%parental.degs.ids)
#also check for deg GO terms 
brass.cultivated.GO.up = topGO_wrapper(geneScores = degs.brass.cultivated,
                                      geneScoresDE = T,
                                      geneScoresDirection = "Up",
                                      GOmapping = GOmapping.brass,
                                      algorithm = "weight01",
                                      statistic = "fisher",
                                      nodeSize = 10,
                                      discretisedDE = T,
                                      p = 0.05)
write.csv(brass.cultivated.GO.up$consolidated_result, 
          file = "Analysis/RNAseq/Tables/brapas_GO_cultivatedbias.csv", row.names = FALSE)
#44 GO terms up inc e.g. responseto bacteria, response to chitin
brass.cultivated.GO.down = topGO_wrapper(geneScores = degs.brass.cultivated,
                                        geneScoresDE = T,
                                        geneScoresDirection = "Down",
                                        GOmapping = GOmapping.brass,
                                        algorithm = "weight01",
                                        statistic = "fisher",
                                        nodeSize = 10,
                                        discretisedDE = T,
                                        p = 0.05)
write.csv(brass.treatment.GO.down$consolidated_result, 
          file = "Analysis/RNAseq/Tables/brapas_GO_wildbias.csv", row.names = FALSE)
#20 down

#results: interaction
degs.brass.interaction = results(dds.gene.deg, 
                                name="treatmentControl.domesticatedCultivated",     
                                alpha = 0.05,
                                lfcThreshold = log2(1))
summary(degs.brass.interaction) #25 degs
degs.brass.interaction.ids = rownames(subset(degs.brass.interaction, padj<=0.05))
degs.brass.interaction.N = length(degs.brass.interaction.ids)
#none of these degs are shared with parental degs
table(degs.brass.interaction.ids%in%parental.degs.ids)
#also check for parental deg GO terms 
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
#16 GO terms

#for the interaction terms, we also plot the output to understand what exactly is going on
#get the 12 terms with lowest adjusted p value
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
  stat_summary(aes(group = domesticated), fun = mean, geom = "path") +
  stat_summary(aes(color = domesticated), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(color = domesticated), fun = mean, geom = "point", size = 4) +
  geom_point(aes(color = domesticated), size = 2) +
  scale_color_manual(values = brewer.pal(6,"Set3")[c(5,6)]) +
  facet_wrap(~gene, scales = "free")

#save plot
ggsave(brass.intplot, 
       filename = "brapas_interaction_deg_plots.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  40, height = 25, units = "cm")

numdegs = c(degs.brass.cultivated.Nup, degs.brass.cultivated.Ndown, 
            degs.brass.treatment.Nup, degs.brass.treatment.Ndown, degs.brass.interaction.N)
numGO = c(nrow(brass.cultivated.GO.up$consolidated_result), nrow(brass.cultivated.GO.down$consolidated_result),
          nrow(brass.treatment.GO.up$consolidated_result), nrow(brass.treatment.GO.down$consolidated_result),
          nrow(brass.interaction.GO$consolidated_result))

brapas.output = data.frame(DEGs=numdegs, GO_terms=numGO, 
                                         row.names = c("Domesticated_bias","Wild_bias","Unstressed_bias","Stressed_bias","Interaction"))

write.csv(brapas.output, file = "Analysis/RNAseq/Tables/brapas_summary.csv")


#### wild vs cultivated brapa interaction norm analysis ####
#one thing we'd like to know: do genes generally gain or lose plasticity after domestication?
#get just cultivated brapa data
metadata.brass.cultivated = subset(metadata.brass.subset, domesticated == "Cultivated")
brass.gene.counts.clean.cultivated = brass.gene.counts.clean[,as.character(metadata.brass.cultivated$sample)]
#get just wild brapa data
metadata.brass.wild = subset(metadata.brass.subset, domesticated == "Wild")
brass.gene.counts.clean.wild = brass.gene.counts.clean[,as.character(metadata.brass.wild$sample)]
#generate fold changes for cultivated
dds.gene.cultivated = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.cultivated,
                                          colData = metadata.brass.cultivated,
                                          design = as.formula(~treatment))
dds.gene.deg.cultivated = DESeq(dds.gene.cultivated, fitType = "parametric", betaPrior = FALSE)
degs.cultivated = results(dds.gene.deg.cultivated,
                       name="treatment_Control_vs_Wheat",
                       alpha = 0.05,
                       lfcThreshold = log2(1))
#generate fold changes for wild
dds.gene.wild = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.wild,
                                               colData = metadata.brass.wild,
                                               design = as.formula(~treatment))
dds.gene.deg.wild = DESeq(dds.gene.wild, fitType = "parametric", betaPrior = FALSE)
degs.wild = results(dds.gene.deg.wild,
                            name="treatment_Control_vs_Wheat",
                            alpha = 0.05,
                            lfcThreshold = log2(1))
#compile foldchange data for the significant interaction terms
foldchanges.subset = data.frame(wild = degs.wild[degs.brass.interaction.ids,"log2FoldChange"],
                                cultivated = degs.cultivated[degs.brass.interaction.ids,"log2FoldChange"], 
                                row.names = degs.brass.interaction.ids)
#cultivated and wild brassica rapa samples do not differ in their plasticity
t.test(x = abs(foldchanges.subset$wild), y = abs(foldchanges.subset$cultivated), paired = TRUE)
wilcox.test(x = abs(foldchanges.subset$wild), y = abs(foldchanges.subset$cultivated), paired = TRUE)
mean(abs(foldchanges.subset$wild), na.rm = T)
mean(abs(foldchanges.subset$cultivated), na.rm = T)
#plot
gg.foldchanges.subset = ggplot(data = reshape2::melt(foldchanges.subset), aes(y = abs(value), x =variable)) +
  geom_boxplot() + 
  #geom_point() +
  labs(x = "group", y = "absolute log2 fold change")
#save plot
ggsave(gg.foldchanges.subset, 
       filename = "brapas_foldchanges_plot.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  30, height = 25, units = "cm")

# another thing we'd like to know: have genes gained or lost plasticity in domestication, 
# and is the direction of that plasticity the same or opposite following domestication?

#pull the relevant data
foldchanges.subset = mutate(foldchanges.subset, 
                            direction=ifelse((sign(raphanistrum)==sign(sativus)),"Equal","Opposite"),
                            magnitude=ifelse((abs(raphanistrum)>abs(sativus)),"Decrease","Increase"))
foldchanges.subset = mutate(foldchanges.subset, 
                            direction=ifelse((sign(wild)==sign(cultivated)),"Equal","Opposite"),
                            magnitude=ifelse((abs(wild)>abs(cultivated)),"Decrease","Increase"))
table.foldchanges = table(select(foldchanges.subset,c("direction","magnitude")))
#distribution of opposite and equal changes is unrelated to magnitudes of changes:
chisq.test(table.foldchanges)
mosaicplot(table.foldchanges)

#more distributed
foldchanges.subset$wildDirection = ifelse(foldchanges.subset$wild>0.5, "up", 
                                          ifelse(foldchanges.subset$wild<(-0.5), "down","neutral"))
foldchanges.subset$cultivatedDirection = ifelse(foldchanges.subset$cultivated>0.5, "up", 
                                                ifelse(foldchanges.subset$cultivated<(-0.5), "down","neutral"))
table.foldchanges.comp = table(select(foldchanges.subset,c("wildDirection","cultivatedDirection")))
#check
chisq.test(table.foldchanges.comp)
mosaicplot(table.foldchanges.comp)

#oppposite direction is significantly more common than same direction
chisq.test(table(select(foldchanges.subset,c("direction"))))
#to confirm, increase in plasticity is significantly more common than decrease
chisq.test(table(select(foldchanges.subset,c("magnitude"))))

# just out of interest, it would be interesting to know whether these overall patterns also hold when looking at 
# the entire set of genes, not just the ones ID'd by DESeq2
#compile foldchange data for all terms
foldchanges.subset.all = data.frame(wild = degs.wild[,"log2FoldChange"],
                                    cultivated = degs.cultivated[,"log2FoldChange"], 
                                    row.names = row.names(degs.wild))
#when looking at all genes (not just ones significant for interaction) cultivated genes are also more plastic overall
t.test(x = abs(foldchanges.subset.all$wild), y = abs(foldchanges.subset.all$cultivated), paired = TRUE)
wilcox.test(x = abs(foldchanges.subset.all$wild), y = abs(foldchanges.subset.all$cultivated), paired = TRUE)
mean(abs(foldchanges.subset.all$wild), na.rm = T)
mean(abs(foldchanges.subset.all$cultivated), na.rm = T)
foldchanges.subset.all = mutate(foldchanges.subset.all, 
                                direction=ifelse((sign(wild)==sign(cultivated)),"Equal","Opposite"),
                                magnitude=ifelse((abs(wild)>abs(cultivated)),"Decrease","Increase"))
#when looking across all genes, neither equal nor reversed directions are more common than one another
chisq.test(table(select(foldchanges.subset.all,c("direction"))))
#however, increase in plasticity is definitely more common even across all genes
chisq.test(table(select(foldchanges.subset.all,c("magnitude"))))







#### wild brapa vs other wilds ####

#now compare wild brapa vs other wild
#subset
metadata.brass.wilds = subset(metadata.brass.clean, domesticated == "Wild")
brass.gene.counts.clean.wilds = brass.gene.counts.clean[,as.character(metadata.brass.wilds$sample)]
#check that metadata and count matrices conform
table(metadata.brass.wilds$sample %in% colnames(brass.gene.counts.clean.wilds))
brass.gene.counts.clean.wilds = brass.gene.counts.clean.wilds[,as.character(metadata.brass.wilds$sample)]
#add column to check whether ancestor of domesticated or not
metadata.brass.wilds$wild.ancestor = (metadata.brass.wilds$species=="Brassica rapa")

#now run model with our variables of interest
dds.gene.wilds = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.wilds,
                                        colData = metadata.brass.wilds,
                                        design = as.formula(~treatment+wild.ancestor+treatment*wild.ancestor))
dds.gene.deg.wilds = DESeq(dds.gene.wilds, fitType = "parametric", betaPrior = FALSE)
# Check for outliers: none are apparent
boxplot(log10(assays(dds.gene.deg.wilds)[["cooks"]]), range=0, las=2)
boxplot(log10(assays(dds.gene.deg.wilds)[["counts"]]), range=0, las=2)

#save for later
save(dds.gene.deg.wilds, file = "Analysis/RNAseq/Tables/brassica_wilds_deseq.R")

resultsNames(dds.gene.deg.wilds)

#results: wheat vs control
degs.brass.wilds.treatment = results(dds.gene.deg.wilds, 
                               name="treatment_Control_vs_Wheat",     
                               alpha = 0.05,
                               lfcThreshold = log2(1))
degs.brass.wilds.treatment.Nup = nrow(subset(degs.brass.wilds.treatment, padj<=0.05 & log2FoldChange>0)) #16 up in control
degs.brass.wilds.treatment.Ndown = nrow(subset(degs.brass.wilds.treatment, padj<=0.05 & log2FoldChange<0)) #22 up in wheat
degs.brass.wilds.treatment.ids = rownames(subset(degs.brass.wilds.treatment, padj<=0.05))
#1 of these degs shared with parental degs
table(degs.brass.wilds.treatment.ids%in%parental.degs.ids)
#also check for  GO terms 
brass.wilds.treatment.GO.up = topGO_wrapper(geneScores = degs.brass.wilds.treatment,
                                      geneScoresDE = T,
                                      geneScoresDirection = "Up",
                                      GOmapping = GOmapping.brass,
                                      algorithm = "weight01",
                                      statistic = "fisher",
                                      nodeSize = 10,
                                      discretisedDE = T,
                                      p = 0.05)
write.csv(brass.wilds.treatment.GO.up$consolidated_result, 
          file = "Analysis/RNAseq/Tables/brapas_wilds_GO_controlbias.csv", row.names = FALSE)
#14 GO terms up
brass.wilds.treatment.GO.down = topGO_wrapper(geneScores = degs.brass.wilds.treatment,
                                        geneScoresDE = T,
                                        geneScoresDirection = "Down",
                                        GOmapping = GOmapping.brass,
                                        algorithm = "weight01",
                                        statistic = "fisher",
                                        nodeSize = 10,
                                        discretisedDE = T,
                                        p = 0.05)
write.csv(brass.wilds.treatment.GO.down$consolidated_result, 
          file = "Analysis/RNAseq/Tables/brapas_wilds_GO_wheatbias.csv", row.names = FALSE)
#20 GO terms down

#results: cultivated vs wild
degs.brass.wilds.cultivated = results(dds.gene.deg.wilds, 
                                name="wild.ancestorTRUE",     
                                alpha = 0.05,
                                lfcThreshold = log2(1))
summary(degs.brass.wilds.cultivated) #~2000 degs
degs.brass.wilds.cultivated.Nup = nrow(subset(degs.brass.wilds.cultivated, padj<=0.05 & log2FoldChange>0)) #1393 up in cultivar
degs.brass.wilds.cultivated.Ndown = nrow(subset(degs.brass.wilds.cultivated, padj<=0.05 & log2FoldChange<0)) #771 up in wild
degs.brass.wilds.cultivated.ids = rownames(subset(degs.brass.wilds.cultivated, padj<=0.05))
#41 of these degs are shared with parental degs
table(degs.brass.wilds.cultivated.ids%in%parental.degs.ids)
#also check for deg GO terms 
brass.wilds.cultivated.GO.up = topGO_wrapper(geneScores = degs.brass.wilds.cultivated,
                                       geneScoresDE = T,
                                       geneScoresDirection = "Up",
                                       GOmapping = GOmapping.brass,
                                       algorithm = "weight01",
                                       statistic = "fisher",
                                       nodeSize = 10,
                                       discretisedDE = T,
                                       p = 0.05)
write.csv(brass.wilds.cultivated.GO.up$consolidated_result, 
          file = "Analysis/RNAseq/Tables/brapas_wilds_GO_cultivatedbias.csv", row.names = FALSE)
#44 GO terms up inc e.g. responseto bacteria, response to chitin
brass.wilds.cultivated.GO.down = topGO_wrapper(geneScores = degs.brass.wilds.cultivated,
                                         geneScoresDE = T,
                                         geneScoresDirection = "Down",
                                         GOmapping = GOmapping.brass,
                                         algorithm = "weight01",
                                         statistic = "fisher",
                                         nodeSize = 10,
                                         discretisedDE = T,
                                         p = 0.05)
write.csv(brass.wilds.treatment.GO.down$consolidated_result, 
          file = "Analysis/RNAseq/Tables/brapas_wilds_GO_wildbias.csv", row.names = FALSE)
#20 down


#results: interaction
degs.brass.wilds.interaction = results(dds.gene.deg.wilds, 
                                      name="treatmentControl.wild.ancestorTRUE",     
                                      alpha = 0.05,
                                      lfcThreshold = log2(1))
summary(degs.brass.wilds.interaction) #73 degs
degs.brass.wilds.interaction.ids = rownames(subset(degs.brass.wilds.interaction, padj<=0.05))
degs.brass.wilds.interaction.N = length(degs.brass.wilds.interaction.ids)
#2 of these degs are shared with parental degs
table(degs.brass.wilds.interaction.ids%in%parental.degs.ids)
#also check for GO terms
brass.wilds.interaction.GO = topGO_wrapper(geneScores = degs.brass.wilds.interaction,
                                          geneScoresDE = T,
                                          geneScoresDirection = NA,
                                          GOmapping = GOmapping.brass,
                                          algorithm = "weight01",
                                          statistic = "fisher",
                                          nodeSize = 10,
                                          discretisedDE = T,
                                          p = 0.05)
nrow(brass.wilds.interaction.GO$consolidated_result)
#33 GO terms

#for the interaction terms, we also plot the output to understand what exactly is going on
#get the 12 terms with lowest adjusted p value
degs.brass.wilds.interaction.signif = subset(degs.brass.wilds.interaction, padj<0.05)
interestgeneswild = row.names(degs.brass.wilds.interaction.signif)[order(degs.brass.wilds.interaction.signif$log2FoldChange,decreasing = T)[1:12]]
#for each gene, extract the counts for plotting and label by gene
for(i in 1:length(interestgeneswild)){
  #if first element, instantiate frame, otherwise rbind to frame
  if(i==1) {
    brass.intplotdata = plotCounts(dds.gene.deg.wilds, gene=interestgeneswild[i], intgroup=c("wild.ancestor","treatment"), returnData = T)
    brass.intplotdata$gene = interestgeneswild[i]
  } else {
    addrows = data.frame(plotCounts(dds.gene.deg.wilds, gene=interestgeneswild[i], intgroup=c("wild.ancestor","treatment"), returnData = T), 
                         gene = interestgeneswild[i])
    brass.intplotdata = rbind(brass.intplotdata, addrows)
  }
}

#plot with ggplot facet wrap
brass.intplot.wilds= ggplot(brass.intplotdata, aes(x = treatment, y = count)) +
  stat_summary(aes(group = wild.ancestor), fun = mean, geom = "path") +
  stat_summary(aes(color = wild.ancestor), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(color = wild.ancestor), fun = mean, geom = "point", size = 4) +
  geom_point(aes(color = wild.ancestor), size = 2) +
  scale_color_manual(values = brewer.pal(6,"Set3")[c(5,6)]) +
  facet_wrap(~gene, scales = "free") +
  labs(color = "Rapa")

#save plot
ggsave(brass.intplot.wilds, 
       filename = "brassica_interaction_deg_plots_wilds.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  40, height = 25, units = "cm")

#data for this need to be instantiated above
wildsnumdegs = c(degs.brass.wilds.cultivated.Nup, degs.brass.wilds.cultivated.Ndown,
                 degs.brass.wilds.treatment.Nup, degs.brass.wilds.treatment.Ndown, degs.brass.wilds.interaction.N)
wildsnumGO = c(nrow(brass.wilds.cultivated.GO.up$consolidated_result), nrow(brass.wilds.cultivated.GO.down$consolidated_result),
               nrow(brass.wilds.treatment.GO.up$consolidated_result), nrow(brass.wilds.treatment.GO.down$consolidated_result),
               nrow(brass.wilds.interaction.GO$consolidated_result))

brass.wilds.output = data.frame(DEGs=wildsnumdegs, GO_terms=wildsnumGO,
                               row.names = c("Domesticated_bias","Wild_bias","Unstressed_bias","Stressed_bias","Interaction"))

write.csv(brass.wilds.output, file = "Analysis/RNAseq/Tables/brass_wilds_summary.csv")

#### wilds interaction norm analysis ####

#first compare treatments for brassica rapa and for other wilds separately

#get fold changes in brapa
#extract norm. counts for each combination
metadata.brass.brapawild = subset(metadata.brass.wilds, species == "Brassica rapa")
metadata.brass.otherwilds = subset(metadata.brass.wilds, species != "Brassica rapa")

brass.gene.counts.clean.otherwilds = brass.gene.counts.clean[,as.character(metadata.brass.otherwilds$sample)]
dds.gene.otherwilds = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.otherwilds,
                                             colData = metadata.brass.otherwilds,
                                             design = as.formula(~treatment))
dds.gene.deg.otherwilds = DESeq(dds.gene.otherwilds, fitType = "parametric", betaPrior = FALSE)
degs.otherwilds = results(dds.gene.deg.otherwilds,
                          name="treatment_Control_vs_Wheat",
                          alpha = 0.05,
                          lfcThreshold = log2(1))

brass.gene.counts.clean.brapawild = brass.gene.counts.clean[,as.character(metadata.brass.brapawild$sample)]
dds.gene.brapawild = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.brapawild,
                                               colData = metadata.brass.brapawild,
                                               design = as.formula(~treatment))
dds.gene.deg.brapawild = DESeq(dds.gene.brapawild, fitType = "parametric", betaPrior = FALSE)
degs.brapawild = results(dds.gene.deg.brapawild,
                            name="treatment_Control_vs_Wheat",
                            alpha = 0.05,
                            lfcThreshold = log2(1))

foldchanges.wilds = data.frame(brapawild = degs.brapawild[degs.brass.wilds.interaction.ids,"log2FoldChange"],
                               others = degs.otherwilds[degs.brass.wilds.interaction.ids,"log2FoldChange"], 
                               row.names = degs.brass.wilds.interaction.ids)

#brapa fold changes are significantly greater than for other wilds
t.test(x = abs(foldchanges.wilds$brapawild), y = abs(foldchanges.wilds$others), paired = TRUE)
wilcox.test(x = abs(foldchanges.wilds$brapawild), y = abs(foldchanges.wilds$others), paired = TRUE)
mean(abs(foldchanges.wilds$brapawild), na.rm = T)
mean(abs(foldchanges.wilds$others), na.rm = T)
#plot
gg.foldchanges.wilds = ggplot(data = reshape2::melt(foldchanges.wilds), aes(y = abs(value), x =variable)) +
  geom_boxplot(range = 0) + 
  #geom_point() +
  labs(x = "group", y = "absolute log2 fold change")
gg.foldchanges.wilds
#save plot
ggsave(gg.foldchanges.wilds, 
       filename = "brapawild_wilds_foldchanges_plot.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  30, height = 25, units = "cm")

#fancier but uglier plot
# ggplot(reshape2::melt(as.matrix(foldchanges.wilds)), aes(y = abs(value))) +
#   geom_boxplot(aes(x = rep(c(-3, 3), each = nrow(foldchanges.wilds)), group = Var2), fill = 'steelblue') +
#   geom_point(aes(x = rep(c(-1, 1), each = nrow(foldchanges.wilds))), size = 5) +
#   geom_line(aes(x = rep(c(-1, 1), each = nrow(foldchanges.wilds)), group = Var1))

#also check whether raphanistrum genes generally have same directionality as other wilds
#pull the relevant data
foldchanges.wilds = mutate(foldchanges.wilds, 
                           direction=ifelse((sign(brapawild)==sign(others)),"Equal","Opposite"),
                           magnitude=ifelse((abs(brapawild)<abs(others)),"Decrease","Increase"))
table.foldchanges = table(select(foldchanges.wilds,c("direction","magnitude")))
#distribution of opposite and equal changes is unrelated to magnitudes of changes:
chisq.test(table.foldchanges)
mosaicplot(table.foldchanges)

#also do more sophisticated comparison of directionality
foldchanges.wilds.all$brapaDirection = ifelse(foldchanges.wilds.all$brapawild>0.5, "up", 
                                              ifelse(foldchanges.wilds.all$brapawild<(-0.5), "down","neutral"))
foldchanges.wilds.all$otherWildsDirection = ifelse(foldchanges.wilds.all$otherwilds>0.5, "up", 
                                                   ifelse(foldchanges.wilds.all$otherwilds<(-0.5), "down","neutral"))
table.foldchanges.wilds.comp = table(select(foldchanges.wilds.all,c("brapaDirection","otherWildsDirection")))
chisq.test(table.foldchanges.wilds.comp)
mosaicplot(table.foldchanges.wilds.comp)

#test differences in plasticity across all genes
t.test(x = abs(foldchanges.wilds$brapawild), y = abs(foldchanges.wilds$others), paired = TRUE)
wilcox.test(x = abs(foldchanges.wilds$brapawild), y = abs(foldchanges.wilds$others), paired = TRUE)

#oppposite direction is significantly more common than same direction, as in raphanus wilds
chisq.test(table(select(foldchanges.wilds,c("direction"))))
#to confirm, brapa wild genes are more likely to exhibit increased plasticity relative to other wilds
chisq.test(table(select(foldchanges.wilds,c("magnitude"))))


# repeat the entire set of genes, not just the ones ID'd by DESeq2
#compile foldchange data for all terms
foldchanges.wilds.all = data.frame(brapawild = degs.brapawild[,"log2FoldChange"],
                                    otherwilds = degs.otherwilds[,"log2FoldChange"], 
                                    row.names = row.names(degs.brapawild))
#when looking at all genes (not just ones significant for interaction) brapa genes are also more plastic overall
t.test(x = abs(foldchanges.wilds.all$brapawild), y = abs(foldchanges.wilds.all$otherwilds), paired = TRUE)
wilcox.test(x = abs(foldchanges.wilds.all$brapawild), y = abs(foldchanges.wilds.all$otherwilds), paired = TRUE)
mean(abs(foldchanges.wilds.all$brapawild), na.rm = T)
mean(abs(foldchanges.wilds.all$otherwilds), na.rm = T)
foldchanges.wilds.all = mutate(foldchanges.wilds.all, 
                                direction=ifelse((sign(brapawild)==sign(otherwilds)),"Equal","Opposite"),
                                magnitude=ifelse((abs(brapawild)>abs(otherwilds)),"Increase","Decrease"))
#when looking across all genes, opposite directionality between brapa and other wilds is more common than chance
chisq.test(table(select(foldchanges.wilds.all,c("direction"))))
#increase in plasticity is definitely more common even across all genes
chisq.test(table(select(foldchanges.wilds.all,c("magnitude"))))

beepr::beep(3)
