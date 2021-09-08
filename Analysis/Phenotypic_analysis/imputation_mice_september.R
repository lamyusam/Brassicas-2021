library("mice")
library("missForest")
library("tidyverse")
library("broom.mixed")

setwd("/home/benjamin/Documents/Brassicas_repo")
#source("Analysis/Phenotypic_analysis/morphology_modeltesting_september.R")

#create object for mice imputation (removing root shoot ratio since this is not a direct measurement)
#phenodata.gen2.clean.brass.test = select(phenodata.gen2.clean.brass, measure.vars.brass) %>% select(-c("Log_root_shoot_ratio"))

#first select predictors- we don't want to impute between variables unless they're correlated (in this case by at least 0.5)
imp = mice(phenodata.gen2.clean.brass, pred=quickpred(phenodata.gen2.clean.brass, mincor=.5), print=F, seed=123)
imp$pred
#some traits don't have a high enough correlation with any other trait and so won't be imputed: 
#Height & Leaf Length 2326, Num Leaves 1820, and SLA & SLDW 1820 (and of course the non-pheno variables)
rowSums(imp$pred)
#For the remaining variables, predictive mean matching will be used, which is appropriate since they're all numeric:
imp$meth
#Check whether imputed values are plausible (blue=original values, red=imputed values)
stripplot(imp, as.formula(paste0(paste(measure.vars.brass,collapse = "+"),"~.imp")),  pch=20, cex=2)
#yes, imputed values look like real values (not surprising given that we're using PMM!)

#now we can pass the imputed values into a regression model
fit = with(imp, lm(Height_1820 ~ Wild_Dom))
#we could extract the results for one iteration of the imputation using:
summary(fit$analyses[[1]])
#instead, let's pool the results of all 5 imputation iterations
pool.fit = pool(fit)
summary(pool.fit)

#okay, now we know that works- let's do the analyses properly
#begin with brassica
ini.brass = mice(phenodata.gen2.clean.brass, pred=quickpred(phenodata.gen2.clean.brass, mincor=.5), print=F, seed=123)
#make sure that secondary variables aren't used for imputation
pred = quickpred(phenodata.gen2.clean.brass, mincor=.5)
pred[ ,row.names(pred)[which(row.names(pred)%in%measure.vars.brass==F)]] = 0
#we can also try using a less conservative but more powerful form of imputation, 
#but my testing suggests this doesn't make much of a difference
# meth = ini.brass$method
# meth[which(meth=="pmm")]="norm"
#refit
imp.brass = mice(phenodata.gen2.clean.brass, pred=pred, method = meth, print=F, seed=123)
#check imputed values still look fine
stripplot(imp.brass, as.formula(paste0(paste(measure.vars.brass,collapse = "+"),"~.imp")),  pch=20, cex=2)

#begin with rapa wild vs dom
imp.brass.subset = filter(imp.brass,Species == "Brassica_rapa")
#loop across all variables and use mice + lmer to model each
for(i in 1:length(measure.vars.brass)){
  
  response.var = measure.vars.brass[i]
  fit = with(imp.brass.subset, lmer(scale(get(response.var))~Environment + Wild_Dom + Environment:Wild_Dom +
                                      (1|Population) + (1|Parental_effects_status)))
  out = summary(pool(fit))
  pvals = data.frame(out$p.value)
  ests = data.frame(out$estimate)
  
  if(i==1){pvals.out = pvals; ests.out = ests}else{
    pvals.out = cbind(pvals.out,pvals); ests.out = cbind(ests.out,ests)}
  
}
#set dimnames
dimnames(pvals.out)=list(c("Intercept","Cultivated_environment","Domesticated_history","Interaction"),measure.vars.brass)
dimnames(ests.out)=list(c("Intercept","Cultivated_environment","Domesticated_history","Interaction"),measure.vars.brass)
#adjust p-vlaues
pvals.out = data.frame(t(apply(pvals.out,1,function(x) p.adjust(x, "BH"))))
#get labels for plot
brassica.subset.coef.labels = paste0(as.matrix(signif(pvals.out,2)),"\n(",as.matrix(signif(ests.out,2)),")")
#coerce df for ggplot and apply FDR correction
brassica.subset.coefs.fdr = reshape2::melt(rownames_to_column(pvals.out)) #%>% mutate(value = p.adjust(value, method = "BH"))
#plot
brassica.subset.traits.fdr.plot = ggplot(data = brassica.subset.coefs.fdr, aes(x = variable, 
                                                                               y = factor(rowname,levels = c("Interaction",
                                                                                                             "Domesticated_history",
                                                                                                             "Cultivated_environment",
                                                                                                             "Intercept")), fill = value)) +
  geom_tile(aes(fill = value),color = "gray", size=.75, width=1, height = 1) +
  geom_text(aes(label=c(brassica.subset.coef.labels), 
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
  ) + coord_fixed()
brassica.subset.traits.fdr.plot
ggsave(brassica.subset.traits.fdr.plot, filename = "phenotypic_lmer_pvals_brassica_subset_imputed.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")



#now repeat for rapa vs other wilds
imp.brass.wilds = filter(imp.brass,Wild_Dom == "Wild")
imp.brass.wilds$data$Progenitor = ifelse(imp.brass.wilds$data$Species=="Brassica_rapa",TRUE,FALSE)
#loop across variables
for(i in 1:length(measure.vars.brass)){
  
  response.var = measure.vars.brass[i]
  fit = with(imp.brass.wilds, lmer(scale(get(response.var))~Environment + Progenitor + Environment:Progenitor +
                                      (1|Population) + (1|Parental_effects_status)))
  out = summary(pool(fit))
  pvals = data.frame(out$p.value)
  ests = data.frame(out$estimate)
  
  if(i==1){pvals.out = pvals; ests.out = ests}else{
    pvals.out = cbind(pvals.out,pvals); ests.out = cbind(ests.out,ests)}
  
}
#set dimnames
dimnames(pvals.out)=list(c("Intercept","Cultivated_environment","Wild_progenitor","Interaction"),measure.vars.brass)
dimnames(ests.out)=list(c("Intercept","Cultivated_environment","Wild_progenitor","Interaction"),measure.vars.brass)
#adjust p-values across rows
pvals.out = data.frame(t(apply(pvals.out,1,function(x) p.adjust(x, "BH"))))
#get labels for plot
brassica.wilds.coef.labels = paste0(as.matrix(signif(pvals.out,2)),"\n(",as.matrix(signif(ests.out,2)),")")
#coerce df for ggplot
brassica.wilds.coefs.fdr = reshape2::melt(rownames_to_column(pvals.out))
#plot
brassica.wilds.traits.fdr.plot = ggplot(data = brassica.wilds.coefs.fdr, 
                                        aes(x = variable, 
                                            y = factor(rowname,levels = c("Interaction",
                                                                          "Wild_progenitor",
                                                                          "Cultivated_environment",
                                                                          "Intercept")), fill = value)) +
  geom_tile(aes(fill = value),color = "gray", size=.75, width=1, height = 1) +
  geom_text(aes(label=c(brassica.wilds.coef.labels), 
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
brassica.wilds.traits.fdr.plot
ggsave(brassica.wilds.traits.fdr.plot, filename = "phenotypic_lmer_pvals_brassica_wilds_imputed.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")


#now for raphanus
ini.raph = mice(phenodata.gen2.clean.raph, pred=quickpred(phenodata.gen2.clean.raph, mincor=.5), print=F, seed=123)
#make sure that secondary variables aren't used for imputation
pred = quickpred(phenodata.gen2.clean.raph, mincor=.5)
pred[ ,row.names(pred)[which(row.names(pred)%in%measure.vars.raph==F)]] = 0
#refit
imp.raph = mice(phenodata.gen2.clean.raph, pred=pred, print=F, seed=123)

#repeat for raphanus sativus vs raphanus raphanistrum
imp.raph.subset = filter(imp.raph, Species %in% c("Raphanus_sativus","Raphanus_raphanistrum"))

#loop over variables
for(i in 1:length(measure.vars.raph)){
  
  response.var = measure.vars.raph[i]
  fit = with(imp.raph.subset, lmer(scale(get(response.var))~Environment + Wild_Dom + Environment:Wild_Dom +
                                      (1|Population) + (1|Parental_effects_status)))
  out = summary(pool(fit))
  pvals = data.frame(out$p.value)
  ests = data.frame(out$estimate)
  
  if(i==1){pvals.out = pvals; ests.out = ests}else{
    pvals.out = cbind(pvals.out,pvals); ests.out = cbind(ests.out,ests)}
  
}
#set dimnames
dimnames(pvals.out)=list(c("Intercept","Cultivated_environment","Domesticated_history","Interaction"),measure.vars.raph)
dimnames(ests.out)=list(c("Intercept","Cultivated_environment","Domesticated_history","Interaction"),measure.vars.raph)
#adjust p-vlaues
pvals.out = data.frame(t(apply(pvals.out,1,function(x) p.adjust(x, "BH"))))
#get labels for plot
raphanus.subset.coef.labels = paste0(as.matrix(signif(pvals.out,2)),"\n(",as.matrix(signif(ests.out,2)),")")
#coerce df for ggplot and apply FDR correction
raphanus.subset.coefs.fdr = reshape2::melt(rownames_to_column(pvals.out)) #%>% mutate(value = p.adjust(value, method = "BH"))
#plot
raphanus.subset.traits.fdr.plot = ggplot(data = raphanus.subset.coefs.fdr, aes(x = variable, 
                                                                               y = factor(rowname,levels = c("Interaction",
                                                                                                             "Domesticated_history",
                                                                                                             "Cultivated_environment",
                                                                                                             "Intercept")),  fill = value)) +
  geom_tile(aes(fill = value),color = "gray", size=.75, width=1, height = 1) +
  geom_text(aes(label=c(raphanus.subset.coef.labels), 
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
raphanus.subset.traits.fdr.plot
ggsave(raphanus.subset.traits.fdr.plot, filename = "phenotypic_lmer_pvals_raphanus_subset_imputed.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")



#now repeat for raphanistrum vs other wilds
imp.raph.wilds = filter(imp.raph, Species!="Raphanus_sativus")
imp.raph.wilds$data$Progenitor = ifelse(imp.raph.wilds$data$Species=="Raphanus_raphanistrum",TRUE,FALSE)
#loop across variables
for(i in 1:length(measure.vars.raph)){
  
  response.var = measure.vars.raph[i]
  fit = with(imp.raph.wilds, lmer(scale(get(response.var))~Environment + Progenitor + Environment:Progenitor +
                                     (1|Population) + (1|Parental_effects_status)))
  out = summary(pool(fit))
  pvals = data.frame(out$p.value)
  ests = data.frame(out$estimate)
  
  if(i==1){pvals.out = pvals; ests.out = ests}else{
    pvals.out = cbind(pvals.out,pvals); ests.out = cbind(ests.out,ests)}
  
}
#set dimnames
dimnames(pvals.out)=list(c("Intercept","Cultivated_environment","Wild_progenitor","Interaction"),measure.vars.raph)
dimnames(ests.out)=list(c("Intercept","Cultivated_environment","Wild_progenitor","Interaction"),measure.vars.raph)
#adjust p-vlaues
pvals.out = data.frame(t(apply(pvals.out,1,function(x) p.adjust(x, "BH"))))
#get labels for plot
raphanus.wilds.coef.labels = paste0(as.matrix(signif(pvals.out,2)),"\n(",as.matrix(signif(ests.out,2)),")")
#coerce df for ggplot and apply FDR correction
raphanus.wilds.coefs.fdr = reshape2::melt(rownames_to_column(pvals.out))
#plot
raphanus.wilds.traits.fdr.plot = ggplot(data = raphanus.wilds.coefs.fdr, 
                                        aes(x = variable, 
                                            y = factor(rowname,levels = c("Interaction",
                                                                          "Wild_progenitor",
                                                                          "Cultivated_environment",
                                                                          "Intercept")), fill = value)) +
  geom_tile(aes(fill = value),color = "gray", size=.75, width=1, height = 1) +
  geom_text(aes(label=c(raphanus.wilds.coef.labels), 
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
raphanus.wilds.traits.fdr.plot
ggsave(raphanus.wilds.traits.fdr.plot, filename = "phenotypic_lmer_pvals_raphanus_wilds_imputed.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")


