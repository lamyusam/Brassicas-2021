library(tidyverse)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(performance)
library(see)
library(DHARMa)
library(reshape2)

setwd("/home/benjamin/Documents/Brassicas_repo")
phenodata.gen2.clean = read.delim("Data/Phenotypic_data/morphology_data_gen2_clean.csv", sep = ",", header = T, row.names = NULL)
colnames(phenodata.gen2.clean)
table(phenodata.gen2.clean$Species)

#refactor to set reference levels
phenodata.gen2.clean$Wild_Dom = factor(phenodata.gen2.clean$Wild_Dom,
                                       c("Wild","Cultivated"))
phenodata.gen2.clean$Environment = factor(phenodata.gen2.clean$Environment,
                                          c("wheat","control"))


# take note of the most relevant measured variables
measure.vars = c("Days_germ",
                 "Height_1820",
                 "Height_2326",
                 "Leaf_length_1820",
                 "Leaf_length_2326",
                 "Num_leaves_1820",
                 "Num_leaves_2326",
                 "Sla_1820",
                 "Sla_2326",
                 "Sldw_1820",
                 "Sldw_2326",
                 "Aboveground_dw",
                 "Root_dw",
                 "Root_to_shoot_ratio"
) 

#narrow to spp
phenodata.gen2.clean.brass = subset(phenodata.gen2.clean, substr(Species,1,3)=="Bra")
phenodata.gen2.clean.raph = subset(phenodata.gen2.clean, substr(Species,1,3)=="Rap")

#Begin by checking distributions for each variable, to check nothing fishy is going on
brassmelt = melt(phenodata.gen2.clean.brass, id.vars = "Population", measure.vars = measure.vars)
gg.brass.hist = ggplot(brassmelt, aes(x=value))+
  geom_histogram() +
  facet_wrap(~variable, scales = "free")
#Repeat for raph
raphmelt = melt(phenodata.gen2.clean.raph, id.vars = "Population", measure.vars = measure.vars)
gg.raph.hist = ggplot(raphmelt, aes(x=value))+
  geom_histogram() +
  facet_wrap(~variable, scales = "free")

#what's going on with root dry weight in these plots? It looks unimodal, but the scale seems off, perhaps indicating outliers
quantile(phenodata.gen2.clean.brass$Root_dw, c(seq(0,0.75,0.25),0.95,1), na.rm=T)
quantile(phenodata.gen2.clean.raph$Root_dw, c(seq(0,0.75,0.25),0.95,1), na.rm=T)
#sure enough, there are definitely outliers here, which we can reasonably attribute to experimental error given the sensitivity of
#the weighing to e.g. residual moisture in the roots. Therefore let's remove these
#quick function to remove values 3+ SDs above the mean
crude_outliercheck = function(x){
  sd = sd(x, na.rm=T)
  mean = mean(x, na.rm=T)
  upper = (mean + sd*3)
  print("Removing values: ");print(paste(x[which(x>upper)]))
  x[which(x>upper)] = NA
  return(x)
}
#apply function to rmeove these obvious outliers
phenodata.gen2.clean.brass$Root_dw = crude_outliercheck(phenodata.gen2.clean.brass$Root_dw)
phenodata.gen2.clean.raph$Root_dw = crude_outliercheck(phenodata.gen2.clean.raph$Root_dw)
#also set root-shoot ratio to NA for these outliers
phenodata.gen2.clean.brass$Root_to_shoot_ratio[which(is.na(phenodata.gen2.clean.brass$Root_dw))] = NA
phenodata.gen2.clean.raph$Root_to_shoot_ratio[which(is.na(phenodata.gen2.clean.raph$Root_dw))] = NA

#Re-check distributions for each variable
brassmelt = melt(phenodata.gen2.clean.brass, id.vars = "Population", measure.vars = measure.vars)
gg.brass.hist = ggplot(brassmelt, aes(x=value))+
  geom_histogram() +
  facet_wrap(~variable, scales = "free")
#Repeat for raph
raphmelt = melt(phenodata.gen2.clean.raph, id.vars = "Population", measure.vars = measure.vars)
gg.raph.hist = ggplot(raphmelt, aes(x=value))+
  geom_histogram() +
  facet_wrap(~variable, scales = "free")
#root dry weight looks better now, but root-shoot ratio is still very skewed- not surprising for a ratio like this
#but if we log the ratio, it looks much better:
hist(log(phenodata.gen2.clean.brass$Root_to_shoot_ratio)); hist(log(phenodata.gen2.clean.raph$Root_to_shoot_ratio))
#okay, let's create a logged ratio
phenodata.gen2.clean.brass$Log_root_shoot_ratio = log(phenodata.gen2.clean.brass$Root_to_shoot_ratio)
phenodata.gen2.clean.raph$Log_root_shoot_ratio = log(phenodata.gen2.clean.raph$Root_to_shoot_ratio)
#update the measure.vars object (also drop traits that remain clearly unimodal)
measure.vars.brass = c(measure.vars[!(measure.vars %in% c("Days_germ","Root_to_shoot_ratio"))],"Log_root_shoot_ratio")
measure.vars.raph = c(measure.vars[!(measure.vars %in% c("Days_germ","Root_to_shoot_ratio"))],"Log_root_shoot_ratio")

#Now we need to try and fit distributions for each... 
#Tom advises not overcomplicating things. For each trait, we'll try a normal distribution and a gamma
for(i in 1:length(measure.vars.brass)){
  
  response_var = measure.vars.brass[i]
  
  phenodata.gen2.clean.brass.omit = subset(phenodata.gen2.clean.brass, is.na(get(response_var))==F)
  model.gauss = lmer(scale(get(response_var))~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population), 
                     data=phenodata.gen2.clean.brass.omit, 
                     control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
  simulation.gauss <- simulateResiduals(fittedModel = model.gauss, plot = F)
  #foo = testDispersion(simulation.gauss)
  residtest.gauss = testResiduals(simulation.gauss)
  #return(plot(simulation.gauss))
  out = data.frame(c(residtest.gauss$uniformity$p.value,
                   residtest.gauss$dispersion$p.value,
                   residtest.gauss$outliers$p.value))
  
  if(i==1){outframe.gauss = out}else{outframe.gauss = cbind(outframe.gauss,out)}
  
}
colnames(outframe.gauss)=measure.vars.brass ; row.names(outframe.gauss) = c("Uniformity","Dispersion","Outliers")
#based on this, the traits mostly look pretty good - even where the p-values are low, the diagnostic plots don't look bad
#let's try again with logged values to be see if that makes anything radically different  
# for(i in 1:length(measure.vars.brass)){
#   
#   response_var = measure.vars.brass[i]
#   
#   phenodata.gen2.clean.brass.omit = subset(phenodata.gen2.clean.brass, is.na(get(response_var))==F)
#   model.log = glmer(get(response_var)~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population), 
#                      data=phenodata.gen2.clean.brass.omit, 
#                      control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
#                      family = gaussian(link="log"))
#   simulation.log <- simulateResiduals(fittedModel = model.log, plot = F)
#   #foo = testDispersion(simulation.log)
#   residtest.log = testResiduals(simulation.log)
#   #return(plot(simulation.log))
#   out = data.frame(c(residtest.log$uniformity$p.value,
#                      residtest.log$dispersion$p.value,
#                      residtest.log$outliers$p.value))
#   
#   if(i==1){outframe.log = out}else{outframe.log = cbind(outframe.log,out)}
#   
# }
# colnames(outframe.log)=measure.vars.brass ; row.names(outframe.log) = c("Uniformity","Dispersion","Outliers")
# outframe.log
#we'd like to try gamma, but doing so frequently throws errors that DHARMa can't handle. Trying an inverse link also throws a lot of errors.
#so for now let's run with simple linear models for brassica, since these seem to work reasonably well. 

#check whether gaussian works for Raphanus traits as well:
for(i in 1:length(measure.vars.raph)){
  
  response_var = measure.vars.raph[i]
  
  phenodata.gen2.clean.raph.omit = subset(phenodata.gen2.clean.raph, is.na(get(response_var))==F)
  model.gauss = lmer(scale(get(response_var))~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population), 
                     data=phenodata.gen2.clean.raph.omit, 
                     control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
  simulation.gauss <- simulateResiduals(fittedModel = model.gauss, plot = F)
  #foo = testDispersion(simulation.gauss)
  residtest.gauss = testResiduals(simulation.gauss)
  #return(plot(simulation.gauss))
  out = data.frame(c(residtest.gauss$uniformity$p.value,
                     residtest.gauss$dispersion$p.value,
                     residtest.gauss$outliers$p.value))
  
  if(i==1){outframe.gauss.raph = out}else{outframe.gauss.raph = cbind(outframe.gauss.raph,out)}
  
}
colnames(outframe.gauss.raph)=measure.vars.raph ; row.names(outframe.gauss.raph) = c("Uniformity","Dispersion","Outliers")
outframe.gauss.raph
#these look good enough
#FORMERLY except for root dry weight which is quite dodgy on the QQplot. 
# phenodata.gen2.clean.raph.omit = subset(phenodata.gen2.clean.raph, is.na(Root_dw)==F)
# hist(scale(phenodata.gen2.clean.raph.omit$Root_dw))
# dw.gamma =  glmer(Root_dw~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population), 
#                    data=phenodata.gen2.clean.raph.omit, 
#                    control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)), family = Gamma())
# simulation.gamma <- simulateResiduals(fittedModel = dw.gamma, plot = F)
# residtest.gamma = testResiduals(simulation.gamma)
# plot(simulation.gamma)
#It also doesn't work with gaussian log or inverse links,nor with gamma, nor poisson (obviously) 

#Number of leaves is a discrete trait, and might be better modeled using a poisson distribution- let's check this
response_var = "Num_leaves_1820"
phenodata.gen2.clean.raph.omit = subset(phenodata.gen2.clean.raph, is.na(get(response_var))==F)
model.pois = glmer(get(response_var)~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population), 
                   data=phenodata.gen2.clean.raph.omit, 
                   control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),family=poisson(link = "log"))
simulation.pois <- simulateResiduals(fittedModel = model.pois, plot = F)
residtest.pois = testResiduals(simulation.pois)
#honestly this doesn't look good, nor does a negative binomial
#perhaps it would be better to drop number of leaves, but we'll simply keep this set of traits Gaussian for now

#okay now it's time to run models for each, remembering to include parental effects as a random factor in all analyses! 
#begin with brassica, comparing wild vs domestic brassica rapa
phenodata.gen2.clean.brass.subset = subset(phenodata.gen2.clean.brass, Species = "Brassica rapa")

for(i in 1:length(measure.vars.brass)){
  
  response_var = measure.vars.brass[i]
  phenodata.gen2.clean.brass.omit = subset(phenodata.gen2.clean.brass.subset, is.na(get(response_var))==F)
  model.gauss = lmer(scale(get(response_var))~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population) + (1|Parental_effects_status), 
                     data=phenodata.gen2.clean.brass.omit, 
                     control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
  out = summary(model.gauss)
  pvals = data.frame(out$coefficients[,"Pr(>|t|)"])
  ests = data.frame(out$coefficients[,"Estimate"])
  
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
ggsave(brassica.subset.traits.fdr.plot, filename = "phenotypic_lmer_pvals_brassica_subset.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")




#now repeat for rapa vs other wilds
#begin with brassica, comparing wild vs domestic brassica rapa
phenodata.gen2.clean.brass.wilds = subset(phenodata.gen2.clean.brass, Wild_Dom == "Wild")
phenodata.gen2.clean.brass.wilds$Progenitor = ifelse(phenodata.gen2.clean.brass.wilds$Species=="Brassica_rapa",TRUE,FALSE)

for(i in 1:length(measure.vars.brass)){
  
  response_var = measure.vars.brass[i]
  phenodata.gen2.clean.brass.omit = subset(phenodata.gen2.clean.brass.wilds, is.na(get(response_var))==F)
  model.gauss = lmer(scale(get(response_var))~Environment + Progenitor + Environment:Progenitor + (1|Population) + (1|Parental_effects_status), 
                     data=phenodata.gen2.clean.brass.omit, 
                     control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
  out = summary(model.gauss)
  pvals = data.frame(out$coefficients[,"Pr(>|t|)"])
  ests = data.frame(out$coefficients[,"Estimate"])
  
  if(i==1){pvals.out = pvals; ests.out = ests}else{
    pvals.out = cbind(pvals.out,pvals); ests.out = cbind(ests.out,ests)}
  
}
#set dimnames
dimnames(pvals.out)=list(c("Intercept","Cultivated_environment","Wild_progenitor","Interaction"),measure.vars.brass)
dimnames(ests.out)=list(c("Intercept","Cultivated_environment","Wild_progenitor","Interaction"),measure.vars.brass)
#adjust p-vlaues
pvals.out = data.frame(t(apply(pvals.out,1,function(x) p.adjust(x, "BH"))))

#get labels for plot
brassica.wilds.coef.labels = paste0(as.matrix(signif(pvals.out,2)),"\n(",as.matrix(signif(ests.out,2)),")")

#coerce df for ggplot and apply FDR correction
brassica.wilds.coefs.fdr = reshape2::melt(rownames_to_column(pvals.out)) #%>% mutate(value = p.adjust(value, method = "BH"))

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
ggsave(brassica.wilds.traits.fdr.plot, filename = "phenotypic_lmer_pvals_brassica_wilds.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")


#repeat for raphanus sativus vs raphanus raphanistrum
phenodata.gen2.clean.raph.subset = subset(phenodata.gen2.clean.raph, Species %in% c("Raphanus_sativus","Raphanus_raphanistrum"))

for(i in 1:length(measure.vars.raph)){
  
  response_var = measure.vars.raph[i]
  phenodata.gen2.clean.raph.omit = subset(phenodata.gen2.clean.raph.subset, is.na(get(response_var))==F)
  model.gauss = lmer(scale(get(response_var))~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population) + (1|Parental_effects_status), 
                     data=phenodata.gen2.clean.raph.omit, 
                     control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
  out = summary(model.gauss)
  pvals = data.frame(out$coefficients[,"Pr(>|t|)"])
  ests = data.frame(out$coefficients[,"Estimate"])
  
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
ggsave(raphanus.subset.traits.fdr.plot, filename = "phenotypic_lmer_pvals_raphanus_subset.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")


#now repeat for raphanistrum vs other wilds
phenodata.gen2.clean.raph.wilds = subset(phenodata.gen2.clean.raph, Species!="Raphanus_sativus")
phenodata.gen2.clean.raph.wilds$Progenitor = ifelse(phenodata.gen2.clean.raph.wilds$Species=="Raphanus_raphanistrum",TRUE,FALSE)

for(i in 1:length(measure.vars.raph)){
  
  response_var = measure.vars.raph[i]
  phenodata.gen2.clean.raph.omit = subset(phenodata.gen2.clean.raph.wilds, is.na(get(response_var))==F)
  model.gauss = lmer(scale(get(response_var))~Environment + Progenitor + Environment:Progenitor + (1|Population) + (1|Parental_effects_status), 
                     data=phenodata.gen2.clean.raph.omit, 
                     control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
  out = summary(model.gauss)
  pvals = data.frame(out$coefficients[,"Pr(>|t|)"])
  ests = data.frame(out$coefficients[,"Estimate"])
  
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
raphanus.wilds.coefs.fdr = reshape2::melt(rownames_to_column(pvals.out)) #%>% mutate(value = p.adjust(value, method = "BH"))

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
ggsave(raphanus.wilds.traits.fdr.plot, filename = "phenotypic_lmer_pvals_raphanus_wilds.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images",
       width =  40, height = 20, units = "cm")
