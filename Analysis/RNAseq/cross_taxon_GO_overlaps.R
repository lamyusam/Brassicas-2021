hyper_overlap = function(x,y,set){
  
  n_A = length(x)
  n_B = length(y)
  n_C = length(set)
  n_A_B = length(intersect(x,y))
  phyper(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = FALSE) 
  
}

setwd("/home/benjamin/Documents/Brassicas_repo/")

load("Data/GO/raph_GOmapping.Rdata")
load("Data/GO/brass_GOmapping.Rdata")

allGO = unique(c(unlist(GOmapping.brass), unlist(GOmapping.raph)))


foo = read.csv("Analysis/RNAseq/Tables/raphanistrum_sativus_GO_wildbias.csv")
bar = read.csv("Analysis/RNAseq/Tables/brapas_GO_wildbias.csv")
hyper_overlap(foo$GO.ID,bar$GO.ID,allGO) #no overlap

foo = read.csv("Analysis/RNAseq/Tables/raphanistrum_sativus_GO_cultivatedbias.csv")
bar = read.csv("Analysis/RNAseq/Tables/brapas_GO_cultivatedbias.csv")
hyper_overlap(foo$GO.ID,bar$GO.ID,allGO) #no overlap

foo = read.csv("Analysis/RNAseq/Tables/raphanistrum_sativus_GO_wheatbias.csv")
bar = read.csv("Analysis/RNAseq/Tables/brapas_GO_wheatbias.csv")
hyper_overlap(foo$GO.ID,bar$GO.ID,allGO)

foo = read.csv("Analysis/RNAseq/Tables/raphanistrum_sativus_GO_controlbias.csv")
bar = read.csv("Analysis/RNAseq/Tables/brapas_GO_controlbias.csv")
hyper_overlap(foo$GO.ID,bar$GO.ID,allGO)



#cultivated vs wild progenitors
foo1 = read.csv("Analysis/RNAseq/Tables/raphanistrum_sativus_GO_wildbias.csv")
bar1 = read.csv("Analysis/RNAseq/Tables/brapas_GO_wildbias.csv")

foo2 = read.csv("Analysis/RNAseq/Tables/raphanistrum_sativus_GO_cultivatedbias.csv")
bar2 = read.csv("Analysis/RNAseq/Tables/brapas_GO_cultivatedbias.csv")

foo3 = unique(c(foo1$GO.ID, foo2$GO.ID))
bar3 = unique(c(bar1$GO.ID, bar2$GO.ID))
hyper_overlap(foo3,bar3,allGO)
intersect(foo3,bar3)
#no overlap
subset(rbind(foo1,foo2), GO.ID %in%intersect(foo3, bar3))$Term


#wild progenitors vs other wilds
foo1 = read.csv("Analysis/RNAseq/Tables/raphanus_wilds_GO_wildbias.csv")
bar1 = read.csv("Analysis/RNAseq/Tables/brapas_wilds_GO_wildbias.csv")

foo2 = read.csv("Analysis/RNAseq/Tables/raphanus_wilds_GO_cultivatedbias.csv")
bar2 = read.csv("Analysis/RNAseq/Tables/brapas_wilds_GO_cultivatedbias.csv")

foo3 = unique(c(foo1$GO.ID, foo2$GO.ID))
bar3 = unique(c(bar1$GO.ID, bar2$GO.ID))
hyper_overlap(foo3,bar3,allGO)
#significant overlap
subset(rbind(foo1,foo2), GO.ID %in%intersect(foo3, bar3))$Term


#stressed vs unstressed (excluding otherwilds)
foo1 = read.csv("Analysis/RNAseq/Tables/raphanistrum_sativus_GO_controlbias.csv")
bar1 = read.csv("Analysis/RNAseq/Tables/brapas_GO_controlbias.csv")

foo2 = read.csv("Analysis/RNAseq/Tables/raphanistrum_sativus_GO_wheatbias.csv")
bar2 = read.csv("Analysis/RNAseq/Tables/brapas_GO_wheatbias.csv")

foo3 = unique(c(foo1$GO.ID, foo2$GO.ID))
bar3 = unique(c(bar1$GO.ID, bar2$GO.ID))
hyper_overlap(foo3,bar3,allGO)
#near-significant overlap
subset(rbind(foo1,foo2), GO.ID %in%intersect(foo3, bar3))$Term



