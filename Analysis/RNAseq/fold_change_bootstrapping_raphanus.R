#here we try repeating the wilds plasticity analysis without combining different species of wild accession together, 
#as this could be supressing the plasticity of the combined samples
#what we're going to do instead is bootstrap three samples of wild rapa vs three samples of another species, multiple times 
#and see what the changes look like there
library(rstatix)


metadata.raph.wilds.combined = subset(metadata.raph, species != "Raphanus sativus")
metadata.raph.wilds.combined$species[which(metadata.raph.wilds.combined$species=="Raphanus raphanistrum mungra")] = "Raphanus raphanistrum munra"

#get fold changes in raphanus
#list all non-rapa wild raphanus species
otherwilds = c("Raphanus sativus var. caudatus","Raphanus raphanistrum munra")
allwilds = c("Raphanus raphanistrum",otherwilds)

#set a log fold-change threshold, if desired
lfc = 1

#pick number of bootstraps to run (keep low for now)
start=Sys.time()
boot = 5
for(sp in 1:length(allwilds)){
  #narrow down to data for focal wild
  focal.wild = as.character(allwilds[sp])
  for(i in 1:boot){
    metadata.raph.focalwild = subset(metadata.raph.wilds.combined, species == focal.wild)
    #pick a random subset of three pairs of wheat and control samples for this species
    pick = sample(which(metadata.raph.focalwild$treatment == "Control"),3)
    #make sure we get matches pairs of samples
    picksamples = stringr::str_sub(metadata.raph.focalwild$label[pick],1,-3)
    metadata.raph.focalwild.sample = subset(metadata.raph.focalwild,substr(label,1,11)%in%picksamples)
    
    #narrow down gene data based on chosen subsample
    raph.gene.counts.clean.focalwild.sample = raph.gene.counts.clean[,as.character(metadata.raph.focalwild.sample$sample)]
    table(metadata.raph.focalwild.sample$sample == colnames(raph.gene.counts.clean.focalwild.sample)) #confirm conformity
    #get DESeq2 output for focal wild
    dds.gene.focalwild = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean.focalwild.sample,
                                                colData = metadata.raph.focalwild.sample,
                                                design = as.formula(~treatment))
    dds.gene.deg.focalwild = DESeq(dds.gene.focalwild, fitType = "parametric", betaPrior = FALSE)
    degs.focalwild = results(dds.gene.deg.focalwild,
                             name="treatment_Control_vs_Wheat",
                             alpha = 0.05, lfcThreshold = lfc)
    
    #get number of degs
    nDEGs = nrow(subset(degs.focalwild,padj<0.05))
    
    #save results of loop
    if(i==1){focalwild.nDEGS = c(nDEGs)}else{
      focalwild.nDEGS = c(focalwild.nDEGS,nDEGs)}
  }
  
  #save data for species
  if(sp == 1){degs.frame.raph = data.frame(focalwild.nDEGS)}else{
    degs.frame.raph = cbind(degs.frame.raph, focalwild.nDEGS)
  }
}
stop=Sys.time()
stop-start
#give colnames to output frame
colnames(degs.frame.raph) = allwilds
beepr::beep(3)
#inconclusive 
#write.csv(degs.frame.raph, file = "Analysis/RNAseq/perwilds_stressDEGtable_raph.csv")
write.csv(degs.frame.raph, file = "Analysis/RNAseq/perwilds_stressDEGtable_raph_lfc1.csv")


#let's try another way and see if the same result comes out: run DESeq2 interaction models for each species vs B rapa
for(i in 1:length(otherwilds)){
  focal = otherwilds[i]
  metadata.raph.focalandraphanistrum = subset(metadata.raph.wilds.combined, species %in% c(focal,"Raphanus raphanistrum"))
  raph.gene.counts.clean.focalandraphanistrum = raph.gene.counts.clean[,as.character(metadata.raph.focalandraphanistrum$sample)]
  #remove spaces so that DESeq doesn't complain
  metadata.raph.focalandraphanistrum$species = str_replace(metadata.raph.focalandraphanistrum$species,pattern = " ",replacement = ".")
  #relevel so that raphica raphanistrum is always the last factor
  metadata.raph.focalandraphanistrum$species = forcats::fct_relevel((metadata.raph.focalandraphanistrum$species),"Raphanus.raphanistrum",after = Inf)
  table(colnames(raph.gene.counts.clean.focalandraphanistrum) == metadata.raph.focalandraphanistrum$sample)
  
  #get DESeq2 output for focal wild
  dds.gene.focalandraphanistrum = DESeqDataSetFromMatrix(countData = raph.gene.counts.clean.focalandraphanistrum,
                                                 colData = metadata.raph.focalandraphanistrum,
                                                 design = as.formula(~treatment+species+treatment*species))
  dds.gene.deg.focalandraphanistrum = DESeq(dds.gene.focalandraphanistrum, fitType = "parametric", betaPrior = FALSE)
  degs.focalandraphanistrum = results(dds.gene.deg.focalandraphanistrum,
                              name="treatmentControl.speciesRaphanus.raphanistrum",
                              alpha = 0.05)
  nUp = nrow(subset(degs.focalandraphanistrum,padj<0.05 & log2FoldChange<0))
  nDown = nrow(subset(degs.focalandraphanistrum,padj<0.05 & log2FoldChange>0))
  if(i==1){interDEGs.raph = data.frame(c(nUp,nDown))}else{
    interDEGs.raph = cbind(interDEGs.raph, data.frame(c(nUp,nDown)))}
}
colnames(interDEGs.raph) = otherwilds; rownames(interDEGs.raph) = c("Up.raphanistrum","Down.raphanistrum")
#this isn't very revealing- there are some meager differences, but nothing conclusive, 
#and nothing that convincingly shows rapa to be more plastic
beepr::beep(3)
write.csv(interDEGs.raph, file = "Analysis/RNAseq/perwilds_interDEGtable_raph.csv")








#melt for purposes of anova
testframe = rownames_to_column(lfcs.frame) %>% reshape2::melt()
res.aov = anova_test(data = testframe, dv = value, wid = rowname, within = variable)
get_anova_table(res.aov)
#pairwise t tests
pwc = testframe %>%
  pairwise_t_test(
    formula = value ~ variable, 
    paired = TRUE,
    ref.group = "Raphanus raphanistrum",
    p.adjust.method = "BH", detailed = T
  )
pwc
#plot
ggplot(testframe, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(x = "Species", y = "LFC unstressed vs stressed") +
  theme_bw() +
  theme(text = element_text(size =12),
        axis.text = element_text(size =12))


#try again with just subset of interaction-significant genes
lfcs.frame.interactiongenes = lfcs.frame[degs.raph.wilds.interaction.ids,]
#melt for purposes of anova
testframe.interactions = rownames_to_column(lfcs.frame.interactiongenes) %>% reshape2::melt()
res.aov = anova_test(data = testframe.interactions, dv = value, wid = rowname, within = variable)
get_anova_table(res.aov)
#pairwise t-tests
pwc = testframe.interactions %>%
  pairwise_t_test(
    formula = value ~ variable, 
    paired = TRUE,
    ref.group = "Raphanus raphanistrum",
    p.adjust.method = "BH",
    detailed = T
  )
pwc
#plot
ggplot(testframe.interactions, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(x = "Species", y = "LFC unstressed vs stressed") +
  theme_bw() +
  theme(text = element_text(size =12),
        axis.text = element_text(size =12))






# #below is to compile data for interaction genes only
# foldchanges.wilds = data.frame(brapawild = degs.brapawild[degs.raph.wilds.interaction.ids,"log2FoldChange"],
#                                focalwild = degs.focalwild[degs.raph.wilds.interaction.ids,"log2FoldChange"],
#                                row.names = degs.raph.wilds.interaction.ids)
#compile fold-change data
foldchanges.wilds = data.frame(brapawild = rowMeans(degsframe.brapawild),
                               focalwild = rowMeans(degsframe.focalwild), 
                               row.names = row.names(degs.brapawild))


#brapa fold changes are significantly greater than for other wilds
t.test(x = abs(foldchanges.wilds$brapawild), y = abs(foldchanges.wilds$focalwild), paired = TRUE)
wilcox.test(x = abs(foldchanges.wilds$brapawild), y = abs(foldchanges.wilds$focalwild), paired = TRUE)
mean(abs(foldchanges.wilds$brapawild), na.rm = T)
mean(abs(foldchanges.wilds$focalwild), na.rm = T)
diffs = data.frame(diff = foldchanges.wilds$brapawild-foldchanges.wilds$focalwild)
#plot
gg.foldchanges.wilds = ggplot(data = diffs, aes(y = diff, x = 1)) +
  geom_boxplot(range = 0) + 
  #geom_point() +
  labs(y = "rapa LFC - cretica LFC")
gg.foldchanges.wilds
