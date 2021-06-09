library(GGally)

data.pairwise.brass = phenodata.gen2.clean.brass[,measure.vars] 

data.pairwise.brass = data.pairwise.brass[,-c(grep("2326",names(data.pairwise.brass)))]

gg.brassica.pairs = ggpairs(data.pairwise.brass, aes(alpha = 0.4), 
                            #upper = list(continuous = "smooth", combo = "box_no_facet", discrete = "facetbar", na ="na")
                            )

ggsave(gg.brassica.pairs, 
       filename = "pairwise_brassicas.png",
       device = "png", path = "/home/benjamin/Documents/Brassicas_repo/Analysis/Phenotypic_analysis/Images/plots_for_Tom/",
       width =  40, height = 40, units = "cm")


data.pairwise.raph = phenodata.gen2.clean.raph[,measure.vars] 

data.pairwise.raph = data.pairwise.raph[,-c(grep("2326",names(data.pairwise.raph)))]

gg.raphanus.pairs = ggpairs(data.pairwise.raph, aes(alpha = 0.4), 
                            #upper = list(continuous = "smooth", combo = "box_no_facet", discrete = "facetbar", na ="na")
)

ggsave(gg.raphanus.pairs, 
       filename = "pairwise_raphanus.png",
       device = "png", path = "/home/benjamin/Documents/Brassicas_repo/Analysis/Phenotypic_analysis/Images/plots_for_Tom/",
       width =  40, height = 40, units = "cm")


data.pairwise.full = phenodata.gen2.clean[,measure.vars]
data.pairwise.full = data.pairwise.full[,-c(grep("2326",names(data.pairwise.full)))]

gg.full.pairs = ggpairs(data.pairwise.full, aes(alpha = 0.4), 
                            #upper = list(continuous = "smooth", combo = "box_no_facet", discrete = "facetbar", na ="na")
)

ggsave(gg.full.pairs, 
       filename = "pairwise_full.png",
       device = "png", path = "/home/benjamin/Documents/Brassicas_repo/Analysis/Phenotypic_analysis/Images/plots_for_Tom/",
       width =  40, height = 40, units = "cm")

