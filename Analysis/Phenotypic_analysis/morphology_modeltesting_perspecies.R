#what we want to do here: use regression models to compare each wild species to the progenitor individually, 
#asking in each case whether there's a significant interaction
#pseudocode:

#list all non-rapa wild brassica species
otherwilds = c("Brassica_cretica","Brassica_incana","Brassica_macrocarpa","Brassica_montana","Brassica_villosa","Brassica_rupestris")
allwilds = c("Brassica_rapa",otherwilds)

for(j in 1:length(otherwilds)){
  focalwild = otherwilds[j]
  
  #subset to just progenitor and chosen wild
  phenodata.gen2.clean.brass.focalpair = subset(phenodata.gen2.clean.brass, Species %in% c("Brassica_rapa",focalwild))
  
  #for each phenotypic variable
  for(i in 1:length(measure.vars.brass)){
    
    response_var = measure.vars.brass[i]
    model.gauss = lmer(scale(get(response_var))~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population), 
                       data=phenodata.gen2.clean.brass.focalpair, 
                       control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
    out = summary(model.gauss)
    #extract info regarding interaction term
    outinter = data.frame(out$coefficients["Environmentcontrol:Wild_DomCultivated",c("Pr(>|t|)","Estimate")]) 
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


