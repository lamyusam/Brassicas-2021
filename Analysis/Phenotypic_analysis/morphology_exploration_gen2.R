library(tidyverse)

setwd("/home/benjamin/Documents/Brassicas_repo")
phenodata.gen2 = read.delim("Data/Phenotypic_data/morphology_data_gen2.txt", sep = ",", header = T)
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
                 "Sp_la", #
                 "Aboveground_dw", #all biomass above soil
                 "Root_dw", #all biomass below soil
                 #"Total_biomass", #sum of above and below ground biomass
                 "Day_flowering", #days from germination to flowering
                 "Root_to_shoot_ratio") #


# convert to numeric where applicable
mutate_at(phenodata.gen2, c(Day_flowering,
                            Num_leaves,
                            Height,
                            sla,
                            sldw,
                            leaf_length,
                            ), )



## change class of some columns
data$Species<-as.character(data$Species)
data$Day_flowering<-as.numeric(as.character(data$Day_flowering))
data$Num_leaves<-as.numeric(as.character(data$Num_leaves))
data$Height<-as.numeric(as.character(data$Height))
data$sla<-as.numeric(as.character(data$sla))
data$sldw<-as.numeric(as.character(data$sldw))
data$leaf_length<-as.numeric(as.character(data$leaf_length))