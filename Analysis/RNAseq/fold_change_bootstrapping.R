#here we try repeating the wilds plasticity analysis without combining different species of wild accession together, 
#as this could be supressing the plasticity of the combined samples
#what we're going to do instead is bootstrap three samples of wild rapa vs three samples of another species, multiple times 
#and see what the changes look like there
library(rstatix)

#get fold changes in brapa
#list all non-rapa wild brassica species
otherwilds = c("Brassica cretica","Brassica incana","Brassica macrocarpa","Brassica montana","Brassica villosa")
allwilds = c("Brassica rapa",otherwilds)
#NB this excludes Rupestris, because that species has only two replicates

#pick number of bootstraps to run (keep low for now)
boot = 2
for(sp in 1:length(allwilds)){
  #narrow down to data for focal wild
  focal.wild = as.character(allwilds[sp])
  for(i in 1:boot){
    metadata.brass.focalwild = subset(metadata.brass.wilds, species == focal.wild)
    #pick a random subset of three pairs of wheat and control samples for this species
    pick = sample(which(metadata.brass.focalwild$treatment == "Control"),3)
    #make sure we get matches pairs of samples
    picksamples = stringr::str_sub(metadata.brass.focalwild$label[pick],1,-3)
    metadata.brass.focalwild.sample = subset(metadata.brass.focalwild,substr(label,1,11)%in%picksamples)

    #narrow down gene data based on chosen subsample
    brass.gene.counts.clean.focalwild.sample = brass.gene.counts.clean[,as.character(metadata.brass.focalwild.sample$sample)]
    table(metadata.brass.focalwild.sample$sample == colnames(brass.gene.counts.clean.focalwild.sample)) #confirm conformity
    #get DESeq2 output for focal wild
    dds.gene.focalwild = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.focalwild.sample,
                                                colData = metadata.brass.focalwild.sample,
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
  if(sp == 1){degs.frame = data.frame(focalwild.nDEGS)}else{
    degs.frame = cbind(degs.frame, focalwild.nDEGS)
  }
}
#give colnames to output frame
colnames(degs.frame) = allwilds
beepr::beep(3)
#now this is really weird... per this, the group that's most plastic by far is Brassica macrocarpa!?

#let's try another way and see if the same result comes out: run DESeq2 interaction models for each species vs B rapa

for(i in 1:length(otherwilds)){
  focal = otherwilds[i]
  metadata.brass.focalandrapa = subset(metadata.brass.wilds, species %in% c(focal,"Brassica rapa"))
  brass.gene.counts.clean.focalandrapa = brass.gene.counts.clean[,as.character(metadata.brass.focalandrapa$sample)]
  #remove spaces so that DESeq doesn't complain
  metadata.brass.focalandrapa$species = str_replace(metadata.brass.focalandrapa$species,pattern = " ",replacement = ".")
  table(colnames(brass.gene.counts.clean.focalandrapa) == metadata.brass.focalandrapa$sample)
  
  #get DESeq2 output for focal wild
  dds.gene.focalandrapa = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.focalandrapa,
                                                 colData = metadata.brass.focalandrapa,
                                                 design = as.formula(~treatment+species+treatment*species))
  dds.gene.deg.focalandrapa = DESeq(dds.gene.focalandrapa, fitType = "parametric", betaPrior = FALSE)
  degs.focalandrapa = results(dds.gene.deg.focalandrapa,
                              name="treatmentControl.speciesBrassica.rapa",
                              alpha = 0.05)
  nUp = nrow(subset(degs.focalandrapa,padj<0.05 & log2FoldChange<0))
  nDown = nrow(subset(degs.focalandrapa,padj<0.05 & log2FoldChange>0))
  if(i==1){interDEGs = data.frame(c(nUp,nDown))}else{
    interDEGs = cbind(interDEGs, data.frame(c(nUp,nDown)))}
}
colnames(interDEGs) = otherwilds







#melt for purposes of anova
testframe = rownames_to_column(lfcs.frame) %>% reshape2::melt()
res.aov = anova_test(data = testframe, dv = value, wid = rowname, within = variable)
get_anova_table(res.aov)
#pairwise t tests
pwc = testframe %>%
  pairwise_t_test(
    formula = value ~ variable, 
    paired = TRUE,
    ref.group = "Brassica rapa",
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
lfcs.frame.interactiongenes = lfcs.frame[degs.brass.wilds.interaction.ids,]
#melt for purposes of anova
testframe.interactions = rownames_to_column(lfcs.frame.interactiongenes) %>% reshape2::melt()
res.aov = anova_test(data = testframe.interactions, dv = value, wid = rowname, within = variable)
get_anova_table(res.aov)
#pairwise t-tests
pwc = testframe.interactions %>%
  pairwise_t_test(
    formula = value ~ variable, 
    paired = TRUE,
    ref.group = "Brassica rapa",
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
# foldchanges.wilds = data.frame(brapawild = degs.brapawild[degs.brass.wilds.interaction.ids,"log2FoldChange"],
#                                focalwild = degs.focalwild[degs.brass.wilds.interaction.ids,"log2FoldChange"],
#                                row.names = degs.brass.wilds.interaction.ids)
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
