library(tidyverse)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(performance)
library(see)

setwd("/home/benjamin/Documents/Brassicas_repo")
phenodata.gen2 = read.delim("Data/Phenotypic_data/morphology_data_gen2.txt", sep = ",", header = T, row.names = NULL)
colnames(phenodata.gen2)
table(phenodata.gen2$Species)

#remove unused columns
phenodata.gen2 = select(phenodata.gen2, -c("Total_la"))

# change the environment labels to something less verbose
phenodata.gen2$Environment =  ifelse(phenodata.gen2$Environment == "wheat competition", "wheat", "control")

# take note of the most relevant measured variables
measure.vars = c("Days_germ", #days from potting to germination
                 "Height", #max distance from soil to plant
                 "Leaf_length", #base to tip of longest non-cotyledon
                 "Num_leaves", #number of leaves
                 "Sla", #single leaf area of oldest non-cotyledon
                 "Sldw", #single lead dry weight of oldest non-cotyledon
                 #"Sp_la", #specific leaf area
                 "Aboveground_dw", #all biomass above soil
                 "Root_dw" #all biomass below soil
                 #"Total_biomass", #sum of above and below ground biomass
                 #"Day_flowering", #days from germination to flowering
                 #"Root_to_shoot_ratio" #ratio of below to above-ground biomass (I think)
                 ) 

# convert to numeric where applicable
phenodata.gen2.clean =  mutate_at(phenodata.gen2, measure.vars, function(x) as.numeric(as.character(x))) %>%
  select(c("Label_front","Species","Population","Wild_Dom","Environment",measure.vars)) 

#the number of complete cases is very small...
dim(phenodata.gen2.clean)
table(complete.cases(phenodata.gen2.clean))
#...but let's see how things look using just complete cases anyway
phenodata.gen2.clean.compcases = phenodata.gen2.clean[complete.cases(phenodata.gen2.clean),]

##pca for brassicas
# subset to brassicas and scale measure vars
pca.counts = phenodata.gen2.clean.compcases %>% 
  subset(substr(Species,1,8)=="Brassica") %>%
  remove_rownames() %>%
  column_to_rownames("Label_front") %>%
  select(measure.vars) #%>% scale()
#create pca object
data.pca = prcomp(t(pca.counts))
#extract PC data
percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
pca.out = list(values = data.frame(data.pca$rotation),
               percent.var = percent.var)
#connect to identifier data
ggpcadata = pca.out$values %>% 
  rownames_to_column(var = "Label_front") %>%
  left_join(phenodata.gen2.clean.compcases,
            by = "Label_front")
#plot
ggplot(ggpcadata, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
  #geom_text(vjust = -1) +
  xlab(paste0("PC",1,": ",signif(pca.out$percent.var[1]*100, 3),"%")) +
  ylab(paste0("PC",2,": ",signif(pca.out$percent.var[2]*100, 3),"%")) +
  theme_bw() +
  scale_color_manual(name = "Treatment",
                     values = brewer.pal(7, "Paired")) +
  # scale_shape_manual(name = "Treatment",
  #                    values = c(8,15:20)) +
  theme(panel.grid = element_line(color = "grey95"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(face = "bold", size =12))

#check distribution of PC assignations
qqnorm(ggpcadata$PC1)
shapiro.test(ggpcadata$PC1)

PC1.model.gaus = lmer(PC1~Environment + Wild_Dom + Environment:Wild_Dom + (1|Wild_Dom/Population), data=ggpcadata)
check_model(PC1.model.gaus)

PC1.model.gamma = glmer(PC1~Environment + Wild_Dom + Environment:Wild_Dom + (1|Wild_Dom/Population), 
                        data=ggpcadata,family="Gamma",glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
check_model(PC1.model.gamma)
anova(PC1.model.gaus, PC1.model.gamma)
summary(PC1.model.gaus)














##pca for radishes
# subset to brassicas and scale measure vars
pca.counts = phenodata.gen2.clean.compcases %>% 
  subset(substr(Species,1,8)=="Raphanus") %>%
  remove_rownames() %>%
  column_to_rownames("Label_front") %>%
  select(measure.vars) %>%
  scale()
#create pca object
data.pca = prcomp(t(pca.counts))
#extract PC data
percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
pca.out = list(values = data.frame(data.pca$rotation),
               percent.var = percent.var)
#connect to identifier data
ggpcadata = pca.out$values %>% 
  rownames_to_column(var = "Label_front") %>%
  left_join(phenodata.gen2.clean.compcases,
            by = "Label_front")
#plot
ggplot(ggpcadata, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
  #geom_text(vjust = -1) +
  xlab(paste0("PC",1,": ",signif(pca.out$percent.var[1]*100, 3),"%")) +
  ylab(paste0("PC",2,": ",signif(pca.out$percent.var[2]*100, 3),"%")) +
  theme_bw() +
  scale_color_manual(name = "Treatment",
                     values = brewer.pal(7, "Paired")) +
  # scale_shape_manual(name = "Treatment",
  #                    values = c(8,15:20)) +
  theme(panel.grid = element_line(color = "grey95"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(face = "bold", size =12))

