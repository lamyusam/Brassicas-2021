#here we try repeating the wilds plasticity analysis without combining different species of wild accession together, 
#as this could be supressing the plasticity of the combined samples
#what we're going to do instead is bootstrap three samples of wild rapa vs three samples of another species, multiple times 
#and see what the changes look like there
library(rstatix)

metadata.raph.wilds.combined = metadata.raph.wilds
metadata.raph.wilds.combined$species[which(metadata.raph.wilds.combined$species=="Raphanus raphanistrum mungra")] = "Raphanus raphanistrum munra"

#get fold changes in brapa
#list all non-rapa wild brassica species
otherwilds = c("Raphanus sativus var. caudatus","Raphanus raphanistrum munra")
allwilds = c("Raphanus raphanistrum",otherwilds)
#NB this excludes Rupestris, because that species has only two replicates

#pick number of bootstraps to run (keep low for now)
boot = 1
for(sp in 1:length(allwilds)){
  #narrow down to data for focal wild
  focal.wild = as.character(allwilds[sp])
  for(i in 1:boot){
    metadata.raph.focalwild = subset(metadata.raph.wilds, species == focal.wild)
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
                             alpha = 0.05)
    
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
#give colnames to output frame
colnames(degs.frame.raph) = allwilds
beepr::beep(3)
#now this is really weird... per this, the group that's most plastic by far is Brassica macrocarpa!?
write.csv(degs.frame.raph, file = "Analysis/RNAseq/perwilds_stressDEGtable_brass.csv")

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
