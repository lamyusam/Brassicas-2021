library(tidyverse)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(performance)
library(see)
library(DHARMa)

setwd("/home/benjamin/Documents/Brassicas_repo")
phenodata.gen2.clean = read.delim("Data/Phenotypic_data/morphology_data_gen2_clean.csv", sep = ",", header = T, row.names = NULL)
colnames(phenodata.gen2.clean)
table(phenodata.gen2.clean$Species)

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
#drop traits that are clearly unimodal:
measure.vars.brass = measure.vars[!(measure.vars %in% c("Days_germ","Root_dw","Root_to_shoot_ratio"))]
measure.vars.raph = measure.vars[!(measure.vars %in% c("Days_germ","Root_to_shoot_ratio"))]

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
for(i in 1:length(measure.vars.brass)){
  
  response_var = measure.vars.brass[i]
  
  phenodata.gen2.clean.brass.omit = subset(phenodata.gen2.clean.brass, is.na(get(response_var))==F)
  model.log = glmer(get(response_var)~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population), 
                     data=phenodata.gen2.clean.brass.omit, 
                     control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
                     family = gaussian(link="log"))
  simulation.log <- simulateResiduals(fittedModel = model.log, plot = F)
  #foo = testDispersion(simulation.log)
  residtest.log = testResiduals(simulation.log)
  #return(plot(simulation.log))
  out = data.frame(c(residtest.log$uniformity$p.value,
                     residtest.log$dispersion$p.value,
                     residtest.log$outliers$p.value))
  
  if(i==1){outframe.log = out}else{outframe.log = cbind(outframe.log,out)}
  
}
colnames(outframe.log)=measure.vars.brass ; row.names(outframe.log) = c("Uniformity","Dispersion","Outliers")
outframe.log
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
#these look mostly fine, except for root dry weight which is quite dodgy on the QQplot. 
phenodata.gen2.clean.raph.omit = subset(phenodata.gen2.clean.raph, is.na(Root_dw)==F)
hist(scale(phenodata.gen2.clean.raph.omit$Root_dw))
dw.gamma =  glmer(Root_dw~Environment + Wild_Dom + Environment:Wild_Dom + (1|Population), 
                   data=phenodata.gen2.clean.raph.omit, 
                   control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)), family = Gamma())
simulation.gamma <- simulateResiduals(fittedModel = dw.gamma, plot = F)
residtest.gamma = testResiduals(simulation.gamma)
plot(simulation.gamma)
#It also doesn't work with gaussian log or inverse links,nor with gamma, nor poisson (obviously) 