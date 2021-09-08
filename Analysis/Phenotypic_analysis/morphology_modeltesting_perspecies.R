#what we want to do here: use regression models to compare each wild species to the progenitor individually, 
#asking in each case whether there's a significant interaction
#pseudocode:

#list all non-rapa wild brassica species
otherwilds = c("Brassica_cretica","Brassica_incana","Brassica_macrocarpa","Brassica_montana","Brassica_villosa","Brassica_rupestris")
allwilds = c("Brassica_rapa",otherwilds)

for(j in 1:length(otherwilds)){
  focalwild = otherwilds[j]
  
  #subset to just progenitor and chosen wild
  phenodata.gen2.clean.brass.focalpair = subset(phenodata.gen2.clean.brass, (Species %in% c("Brassica_rapa",focalwild)) & Wild_Dom=="Wild")
  phenodata.gen2.clean.brass.focalpair$Progenitor = ifelse(phenodata.gen2.clean.brass.focalpair$Species == "Brassica_rapa",TRUE,FALSE)
  
  #for each phenotypic variable
  for(i in 1:length(measure.vars.brass)){
    
    response_var = measure.vars.brass[i]
    model.gauss = lmer(scale(get(response_var))~Environment + Progenitor + Environment:Progenitor + (1|Population), 
                       data=phenodata.gen2.clean.brass.focalpair, 
                       control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
    out = summary(model.gauss)
    #extract info regarding interaction term
    outinter = data.frame(out$coefficients["Environmentcontrol:ProgenitorTRUE",c("Pr(>|t|)","Estimate")]) 
    #save loop output
    if(i==1){outframe = outinter}else{outframe = cbind(outframe,outinter)}
    
  }
  
  #collate info per spp loop
  dimnames(outframe)=list(c("pval","estimate"),measure.vars.brass)
  if(j==1){pvals.out = outframe["pval",]; ests.out = outframe["estimate",]}else{
    pvals.out = rbind(pvals.out,outframe["pval",]); ests.out = rbind(ests.out,outframe["estimate",])}
}
  
#rename
row.names(pvals.out) = otherwilds ; row.names(ests.out) = otherwilds 
#adjust p-vlaues
pvals.out = data.frame(apply(pvals.out,2,function(x) p.adjust(x, "BH")))

#get labels for plot
brass.wilds.coefs.labels = paste0(as.matrix(signif(pvals.out,2)),"\n(",as.matrix(signif(ests.out,2)),")")

#coerce df for ggplot 
brass.wilds.coefs.fdr = reshape2::melt(rownames_to_column(pvals.out)) 

#plot
brass.wilds.speciesinteractions.fdr.plot = ggplot(data = brass.wilds.coefs.fdr, aes(x = variable,  y = rowname,  fill = value)) +
  geom_tile(aes(fill = value),color = "gray", size=.75, width=1, height = 1) +
  geom_text(aes(label=c(brass.wilds.coefs.labels), 
                lineheight = 0.75, size = 2), show.legend = FALSE) +
  scale_fill_gradientn(colours = colorRampPalette(rev(c("#FFFFFF",brewer.pal(n = 9, name = "Reds")[1:5])),bias=6)(20),
                       breaks = c(0.0,0.05,1),
                       expand = c(0,0),
                       limits = c(0,1),
                       guide = guide_colourbar(barheight = 25,
                                               #title = "p-value\n(adjusted)",
                                               title = "p-value\n(adjusted)",
                                               title.vjust = 2,
                                               frame.colour = "black", 
                                               frame.linewidth = 1.5)) +
  theme_minimal() + 
  theme(#aspect.ratio = 1,
    panel.grid = element_line(size = 0.2, colour = "gray80"),
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank()#,
    #legend.text = element_text(size=12),
    #legend.title = element_text(size=15, face = "bold")
  ) +
  coord_fixed()

brass.wilds.speciesinteractions.fdr.plot
ggsave(brass.wilds.speciesinteractions.fdr.plot, filename = "phenotypic_lmer_pvals_brass_wilds_perspecies.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")




#repeat for raphanus
otherwilds = c("Raphanus_sativus_var_caudatus","Raphanus_raphanistrum_munra")
allwilds = c("Raphanus_raphanistrum",otherwilds)

for(j in 1:length(otherwilds)){
  focalwild = otherwilds[j]
  
  #subset to just progenitor and chosen wild
  phenodata.gen2.clean.raph.focalpair = subset(phenodata.gen2.clean.raph, Species %in% c("Raphanus_raphanistrum",focalwild))
  phenodata.gen2.clean.raph.focalpair$Progenitor = ifelse(phenodata.gen2.clean.raph.focalpair$Species == "Raphanus_raphanistrum",TRUE,FALSE)
  
  #for each phenotypic variable
  for(i in 1:length(measure.vars.raph)){
    
    response_var = measure.vars.raph[i]
    model.gauss = lmer(scale(get(response_var))~Environment + Progenitor + Environment:Progenitor + (1|Population), 
                       data=phenodata.gen2.clean.raph.focalpair, 
                       control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
    out = summary(model.gauss)
    #extract info regarding interaction term
    outinter = data.frame(out$coefficients["Environmentcontrol:ProgenitorTRUE",c("Pr(>|t|)","Estimate")]) 
    #save loop output
    if(i==1){outframe = outinter}else{outframe = cbind(outframe,outinter)}
    
  }
  
  #collate info per spp loop
  dimnames(outframe)=list(c("pval","estimate"),measure.vars.raph)
  if(j==1){pvals.out = outframe["pval",]; ests.out = outframe["estimate",]}else{
    pvals.out = rbind(pvals.out,outframe["pval",]); ests.out = rbind(ests.out,outframe["estimate",])}
}

#rename
row.names(pvals.out) = otherwilds ; row.names(ests.out) = otherwilds 
#adjust p-vlaues
pvals.out = data.frame(apply(pvals.out,2,function(x) p.adjust(x, "BH")))

#get labels for plot
raph.wilds.coefs.labels = paste0(as.matrix(signif(pvals.out,2)),"\n(",as.matrix(signif(ests.out,2)),")")

#coerce df for ggplot 
raph.wilds.coefs.fdr = reshape2::melt(rownames_to_column(pvals.out)) 

#plot
raph.wilds.speciesinteractions.fdr.plot = ggplot(data = raph.wilds.coefs.fdr, aes(x = variable,  y = rowname,  fill = value)) +
  geom_tile(aes(fill = value),color = "gray", size=.75, width=1, height = 1) +
  geom_text(aes(label=c(raph.wilds.coefs.labels), 
                lineheight = 0.75, size = 2), show.legend = FALSE) +
  scale_fill_gradientn(colours = colorRampPalette(rev(c("#FFFFFF",brewer.pal(n = 9, name = "Reds")[1:5])),bias=6)(20),
                       breaks = c(0.0,0.05,1),
                       expand = c(0,0),
                       limits = c(0,1),
                       guide = guide_colourbar(barheight = 25,
                                               #title = "p-value\n(adjusted)",
                                               title = "p-value\n(adjusted)",
                                               title.vjust = 2,
                                               frame.colour = "black", 
                                               frame.linewidth = 1.5)) +
  theme_minimal() + 
  theme(#aspect.ratio = 1,
    panel.grid = element_line(size = 0.2, colour = "gray80"),
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank()#,
    #legend.text = element_text(size=12),
    #legend.title = element_text(size=15, face = "bold")
  ) +
  coord_fixed()

raph.wilds.speciesinteractions.fdr.plot
ggsave(raph.wilds.speciesinteractions.fdr.plot, filename = "phenotypic_lmer_pvals_raph_wilds_perspecies.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")

