#here we try repeating the wilds plasticity analysis without combining different species of wild accession together, 
#as this could be supressing the plasticity of the combined samples
#what we're going to do instead is bootstrap three samples of wild rapa vs three samples of another species, multiple times 
#and see what the changes look like there


#get fold changes in brapa
#list all non-rapa wild brassica species
otherwilds = c("Brassica cretica","Brassica incana","Brassica macrocarpa","Brassica montana","Brassica villosa")
allwilds = c("Brassica rapa",otherwilds)
#NB this excludes Rupestris, because that species has only two replicates

#pick number of bootstraps to run (keep low for now)
boot = 3
#narrow down to data for focal wild
focal.wild = "Brassica macrocarpa"
for(i in 1:boot){
  metadata.brass.focalwild = subset(metadata.brass.wilds, species == focal.wild)
  #pick a random subset of three wheat and three control samples for this species
  pick = c(sample(which(metadata.brass.focalwild$treatment == "Control"),3),
           sample(which(metadata.brass.focalwild$treatment == "Wheat"),3))
  metadata.brass.focalwild.sample = metadata.brass.focalwild[pick,]
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
  #save results of loop for wild
  if(i==1){degsframe.focalwild = data.frame(degs.focalwild$log2FoldChange)}else{
    degsframe.focalwild = cbind(degsframe.focalwild,degs.focalwild$log2FoldChange)}
  
  #repeat process for wild brassica rapa
  #pick a random subset of three wheat and three control samples for this species
  pick = c(sample(which(metadata.brass.brapawild$treatment == "Control"),3),
           sample(which(metadata.brass.brapawild$treatment == "Wheat"),3))
  metadata.brass.brapawild.sample = metadata.brass.brapawild[pick,]
  #narrow down gene data based on chosen subsample
  brass.gene.counts.clean.brapawild.sample = brass.gene.counts.clean[,as.character(metadata.brass.brapawild.sample$sample)]
  table(metadata.brass.brapawild.sample$sample == colnames(brass.gene.counts.clean.brapawild.sample)) #confirm conformity
  #get DESeq2 output for focal wild
  dds.gene.brapawild = DESeqDataSetFromMatrix(countData = brass.gene.counts.clean.brapawild.sample,
                                              colData = metadata.brass.brapawild.sample,
                                              design = as.formula(~treatment))
  dds.gene.deg.brapawild = DESeq(dds.gene.brapawild, fitType = "parametric", betaPrior = FALSE)
  degs.brapawild = results(dds.gene.deg.brapawild,
                           name="treatment_Control_vs_Wheat",
                           alpha = 0.05)
  #save results of loop for wild
  if(i==1){degsframe.brapawild = data.frame(degs.brapawild$log2FoldChange)}else{
    degsframe.brapawild = cbind(degsframe.brapawild,degs.brapawild$log2FoldChange)}
  
}
beepr::beep(3)

#below is to compile data for interaction genes only
foldchanges.wilds = data.frame(brapawild = degs.brapawild[degs.brass.wilds.interaction.ids,"log2FoldChange"],
                               focalwild = degs.focalwild[degs.brass.wilds.interaction.ids,"log2FoldChange"],
                               row.names = degs.brass.wilds.interaction.ids)
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
