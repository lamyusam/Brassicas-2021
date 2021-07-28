library("mice")
library("missForest")
library("tidyverse")
library("sinkr")
library("pcaMethods")
library("factoextra")


# narrow down to non-biomass data and remove 1820 data since these are mostly mutually exclusive with non-1820 data
# also remove days germ, because the distribution is so heavily unimodal
phenodata.gen2.raph.pcadata = phenodata.gen2.clean.raph %>% 
  remove_rownames() %>%
  column_to_rownames(var = "Label_front") %>%
  select(-c("Aboveground_dw","Root_dw","Root_to_shoot_ratio","Days_germ")) %>%  
  select(-grep(names(phenodata.gen2.clean.raph), pattern = "1820",value = T))

# first generate PCA using just complete cases
phenodata.gen2.raph.pcadata.full = phenodata.gen2.raph.pcadata[complete.cases(phenodata.gen2.raph.pcadata),]
#create pca object
data.pca = prcomp(phenodata.gen2.raph.pcadata.full,scale.=T)
#display variable loadings
pca.var.plot.raph = fviz_pca_var(data.pca,
                                 col.var = "contrib", # Color by contributions to the PC
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE ,    # Avoid text overlapping
                                 title = "Raphanus")
#display samples
pca.ind.plot.raph = fviz_pca_ind(data.pca,
                                 label = "none", # hide individual labels
                                 habillage = phenodata.gen2.raph.pcadata.full$Species, # color by species
                                 #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 addEllipses = TRUE) # Concentration ellipses


#then impute data three ways, and for each re-generate PCA to compare to original