library("mice")
library("missForest")
library("tidyverse")
library("sinkr")
library("pcaMethods")
library("factoextra")
library("ggpubr")
library("DHARMa")

setwd("/home/benjamin/Documents/Brassicas_repo")

# narrow down to non-biomass data and remove 1820 data since these are mostly mutually exclusive with non-1820 data
# also remove days germ, because the distribution is so heavily unimodal
phenodata.gen2.raph.pcadata = phenodata.gen2.clean.raph %>% 
  remove_rownames() %>%
  column_to_rownames(var = "Label_front") %>%
  select(measure.vars) %>%
  select(-c("Aboveground_dw","Root_dw","Root_to_shoot_ratio","Days_germ")) %>%  
  select(-grep(names(phenodata.gen2.clean.raph), pattern = "1820",value = T)) 

#also remove any rows that are comprised only of NAs
phenodata.gen2.raph.pcadata = phenodata.gen2.raph.pcadata[rowSums(is.na(phenodata.gen2.raph.pcadata))!=ncol(phenodata.gen2.raph.pcadata),]

# first generate PCA using just complete cases
set.seed(123)
phenodata.gen2.raph.pcadata.full = phenodata.gen2.raph.pcadata[complete.cases(phenodata.gen2.raph.pcadata),]
#create pca object
data.pca = prcomp(phenodata.gen2.raph.pcadata.full,scale.=T)
#display variable loadings
pca.var.plot.raph = fviz_pca_var(data.pca,
                                 col.var = "contrib", # Color by contributions to the PC
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE ,    # Avoid text overlapping
                                 title = "No imputation")
# #display samples
# pca.ind.plot.raph = fviz_pca_ind(data.pca,
#                                  label = "none", # hide individual labels
#                                  habillage = , # color by species
#                                  #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#                                  addEllipses = TRUE) # Concentration ellipses


# now generate PCA using imputed data from mice
set.seed(123)
init = mice(phenodata.gen2.raph.pcadata, m=5, method = 'pmm')
meth = init$method
predM = init$predictorMatrix
imputed.mice = mice(phenodata.gen2.raph.pcadata, method=meth, predictorMatrix=predM, m=5)
imputed.mice = complete(imputed.mice)
#create pca object
data.pca.mice = prcomp(imputed.mice,scale.=T)
data.pca.mice$rotation[,"PC2"] = data.pca.mice$rotation[,"PC2"]*(-1) #we flip PC2 here for consistency with the other plots
#display variable loadings
pca.var.plot.raph.mice = fviz_pca_var(data.pca.mice,
                                 col.var = "contrib", # Color by contributions to the PC
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE,
                                 title = "MICE")
# #display samples
# pca.ind.plot.raph = fviz_pca_ind(data.pca,
#                                  label = "none", # hide individual labels
#                                  habillage = phenodata.gen2.raph.pcadata.full$Species, # color by species
#                                  #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#                                  addEllipses = TRUE) # Concentration ellipses

# now generate PCA using imputed data from dineof
set.seed(123)
imputed.dineof = dineof(Xo = as.matrix(phenodata.gen2.raph.pcadata), delta.rms = 1e-04)
#create pca object
data.pca.dineof = prcomp(imputed.dineof$Xa,scale.=T) 
data.pca.dineof$rotation[,"PC2"] = data.pca.dineof$rotation[,"PC2"]*(-1) #we flip PC2 here for consistency with the other plots
#display variable loadings
pca.var.plot.raph.dineof = fviz_pca_var(data.pca.dineof,
                                      col.var = "contrib", # Color by contributions to the PC
                                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                      repel = TRUE,
                                      title = "DINEOF")
# #display samples
# pca.ind.plot.raph = fviz_pca_ind(data.pca,
#                                  label = "none", # hide individual labels
#                                  habillage = phenodata.gen2.raph.pcadata.full$Species, # color by species
#                                  #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#                                  addEllipses = TRUE) # Concentration ellipses


# now generate PCA using imputed data from pcamethods
set.seed(123)
imputed.methods.pca = pca(phenodata.gen2.raph.pcadata[rowSums(is.na(phenodata.gen2.raph.pcadata))!=
                                                        ncol(phenodata.gen2.raph.pcadata),], nPcs=2, method="ppca")
imputed.methods = completeObs(imputed.methods.pca)
#create pca object
data.pca.methods = prcomp(imputed.methods,scale.=T)
data.pca.methods$rotation[,"PC2"] = data.pca.methods$rotation[,"PC2"]*(-1) #we flip PC2 here for consistency with the other plots
#display variable loadings
pca.var.plot.raph.methods = fviz_pca_var(data.pca.methods,
                                        col.var = "contrib", # Color by contributions to the PC
                                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                        repel = TRUE,
                                        title = "PCAmethods")
# #display samples
# pca.ind.plot.raph = fviz_pca_ind(data.pca,
#                                  label = "none", # hide individual labels
#                                  habillage = phenodata.gen2.raph.pcadata.full$Species, # color by species
#                                  #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#                                  addEllipses = TRUE) # Concentration ellipses

#arrange and save
pca.arrange = ggarrange(pca.var.plot.raph, pca.var.plot.raph.mice, pca.var.plot.raph.dineof, pca.var.plot.raph.methods, common.legend = T)
ggsave(pca.arrange, 
       filename = "pca_imputation_plots.pdf",
       device = "pdf", path = "Analysis/Phenotypic_analysis/Images/",
       width =  40, height = 35, units = "cm")


# now we want to see if using the imputed data allow us to identify a PC-level response to stress/wildness
attachdata = select(phenodata.gen2.clean.raph, c("Label_front","Wild_Dom","Environment","Species","Population"))

#begin with un-imputed data
compframe = data.frame(data.pca$x[,c("PC1","PC2")]) %>% rownames_to_column(var = "Label_front")
compframe = left_join(compframe, attachdata, by = "Label_front")
#for PC1
model.pc1 = lmer(formula = "PC1~Wild_Dom+Environment+(1|Population)", data = compframe)
summary(model.pc1)$coefficients
simulationOutput = simulateResiduals(fittedModel = model.pc1)
plot(simulationOutput)
write.csv(signif(summary(model.pc1)$coefficients, digits = 4), 
          file = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Tables/model.pc1.csv")
#for PC2
model.pc2 = lmer(formula = "PC2~Wild_Dom+Environment+(1|Population)", data = compframe)
summary(model.pc2)$coefficients
simulationOutput = simulateResiduals(fittedModel = model.pc2)
plot(simulationOutput)
write.csv(signif(summary(model.pc2)$coefficients, digits = 4), 
          file = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Tables/model.pc2.csv")


#now with MICE data
miceframe = data.frame(data.pca.mice$x[,c("PC1","PC2")]) %>% rownames_to_column(var = "Label_front")
miceframe = left_join(miceframe, attachdata, by = "Label_front")
#for PC1
model.mice.pc1 = lmer(formula = "PC1~Wild_Dom+Environment+(1|Population)", data = miceframe)
summary(model.mice.pc1)$coefficients
simulationOutput = simulateResiduals(fittedModel = model.mice.pc1)
plot(simulationOutput)
write.csv(signif(summary(model.mice.pc1)$coefficients, digits = 4), 
          file = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Tables/model.mice.pc1.csv")
#for PC2
model.mice.pc2 = lmer(formula = "PC2~Wild_Dom+Environment+(1|Population)", data = miceframe)
summary(model.mice.pc2)$coefficients
simulationOutput = simulateResiduals(fittedModel = model.mice.pc2)
plot(simulationOutput)
write.csv(signif(summary(model.mice.pc2)$coefficients, digits = 4), 
          file = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Tables/model.mice.pc2.csv")


#now with pcamethods data
methodsframe = data.frame(data.pca.methods$x[,c("PC1","PC2")]) %>% rownames_to_column(var = "Label_front")
methodsframe = left_join(methodsframe, attachdata, by = "Label_front")
#for PC1
model.methods.pc1 = lmer(formula = "PC1~Wild_Dom+Environment+(1|Population)", data = methodsframe)
summary(model.methods.pc1)$coefficients
simulationOutput = simulateResiduals(fittedModel = model.methods.pc1)
plot(simulationOutput)
write.csv(signif(summary(model.methods.pc1)$coefficients, digits = 4), 
          file = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Tables/model.methods.pc1.csv")
#for PC2
model.methods.pc2 = lmer(formula = "PC2~Wild_Dom+Environment+(1|Population)", data = methodsframe)
summary(model.methods.pc2)$coefficients
simulationOutput = simulateResiduals(fittedModel = model.methods.pc2)
plot(simulationOutput)
write.csv(signif(summary(model.methods.pc2)$coefficients, digits = 4), 
          file = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Tables/model.methods.pc2.csv")



#also with dineof data
dineofframe = data.frame(data.pca.dineof$x[,c("PC1","PC2")]) %>% rownames_to_column(var = "Label_front")
dineofframe = left_join(dineofframe, attachdata, by = "Label_front")
#for PC1
model.dineof.pc1 = lmer(formula = "PC1~Wild_Dom+Environment+(1|Population)", data = dineofframe)
summary(model.dineof.pc1)$coefficients
simulationOutput = simulateResiduals(fittedModel = model.dineof.pc1)
plot(simulationOutput)
write.csv(signif(summary(model.dineof.pc1)$coefficients, digits = 4), 
          file = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Tables/model.dineof.pc1.csv")
#for PC2
model.dineof.pc2 = lmer(formula = "PC2~Wild_Dom+Environment+(1|Population)", data = dineofframe)
summary(model.dineof.pc2)$coefficients
simulationOutput = simulateResiduals(fittedModel = model.dineof.pc2)
plot(simulationOutput)
write.csv(signif(summary(model.dineof.pc2)$coefficients, digits = 4), 
          file = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Tables/model.dineof.pc2.csv")
