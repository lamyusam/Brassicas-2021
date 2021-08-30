library(RRHO)
library(stringr)
library(tidyverse)
library(pheatmap)
library(DESeq2)

#set working directory
setwd("/home/benjamin/Documents/Brassicas_repo")

#get Blast hits
brass_raph_RBH =  read.delim("Data/RNAseq/raph_brass_RBH.txt", sep = " ", header = F, col.names = c("brass","raph"))
#rename to conform to count matrices
brass_raph_RBH = apply(brass_raph_RBH, MARGIN = 2, FUN = function(x) paste0("gene-",x))

#we want to retain only the hits for which both genes are present in the respective raphanus and brassica count matrices
allnames = c(row.names(brass.gene.counts.clean), row.names(raph.gene.counts.clean))
brass_raph_RBH = data.frame(brass_raph_RBH[rowSums(apply(brass_raph_RBH, MARGIN = 2, FUN = function(x) x %in% allnames))==2,])
#this leaves a good number of genes- 24315

#let's create an arbitrary label for each gene pairing to avoid confusion
brass_raph_RBH$label = paste0("G",1:nrow(brass_raph_RBH))

#load dds objects for raphanus and brassica
brassrapa.dds = load("Analysis/RNAseq/Tables/brassica_rapa_deseq.R") ; brassrapa.dds = get(brassrapa.dds)
brasswilds.dds = load("Analysis/RNAseq/Tables/brassica_wilds_deseq.R") ; brasswilds.dds = get(brasswilds.dds)
raphsubset.dds = load("Analysis/RNAseq/Tables/raphanus_sativus_raphanistrum_deseq.R") ; raphsubset.dds = get(raphsubset.dds)
raphwilds.dds = load("Analysis/RNAseq/Tables/raphanus_wilds_deseq.R") ; raphwilds.dds = get(raphwilds.dds)

#get objects for a comparison of interest, for example wild progenitor vs domesticate
degs.brass.cultivated = results(brassrapa.dds, name="domesticated_Cultivated_vs_Wild", alpha = 0.05)     
degs.raph.cultivated = results(raphsubset.dds, name="domesticated_Cultivated_vs_Wild", alpha = 0.05)   
#for stress in wild progenitor and domesticate
degs.brass.cultivated.stress = results(brassrapa.dds, name="treatment_Control_vs_Wheat", alpha = 0.05)
degs.raph.cultivated.stress = results(raphsubset.dds, name="treatment_Control_vs_Wheat", alpha = 0.05)
#for wild progenitor and other wilds
degs.brass.wilds = results(brasswilds.dds, name="wild.ancestorTRUE", alpha = 0.05)     
degs.raph.wilds = results(raphwilds.dds, name="wild.ancestorTRUE", alpha = 0.05) 
#for stress in wild progenitor and other wilds
degs.brass.wilds.stress = results(brasswilds.dds, name="treatment_Control_vs_Wheat", alpha = 0.05)     
degs.raph.wilds.stress = results(raphwilds.dds, name="treatment_Control_vs_Wheat", alpha = 0.05) 

prep4RRHO = function(brassdegs, raphdegs){

  #remove rows that aren't in the matching table
  brassdegs = subset(brassdegs, row.names(brassdegs) %in% brass_raph_RBH$brass)
  #relabel each output's rownames to match the labels in the Blast output table
  row.names(brassdegs) = brass_raph_RBH$label[match(brass_raph_RBH$brass, row.names(brassdegs))]
  
  #remove rows that aren't in the matching table
  raphdegs = subset(raphdegs, row.names(raphdegs) %in% brass_raph_RBH$raph)
  #relabel each output's rownames to match the labels in the Blast output table
  row.names(raphdegs) = brass_raph_RBH$label[match(brass_raph_RBH$raph, row.names(raphdegs))]
  
  #reorder the frames to match
  raphdegs = raphdegs[row.names(brassdegs),]
  
  #now we should be able to input to RRHO
  #before transforming for RRHO, we might wish to remove the set of genes that have p values close to 1 in both comparisons
  #this may improve the analysis because genes with very low significance are unlikely to rank meaningfully
  #combine two pvalue lists
  both = cbind(brassdegs$pvalue, raphdegs$pvalue) 
  #remove rows that have p-value above threshold in both comparisons
  keep =  row.names(brassdegs)[(rowSums(both > 0.8 | is.na(both))!=2)]
  print(paste0("Removing ",nrow(brassdegs)-length(keep)," genes with low significance in both comparisons."))
  
  #log transformation of p-values
  DEG2RRHO = function(contrast){
    
    RRHO = data.frame(gene = row.names(contrast),
                      value = -log10(contrast$pvalue)*sign(contrast$log2FoldChange))
    #value = contrast$log2FoldChange)
    return(RRHO)
    
  }
  
  brassdegs.RRHO = DEG2RRHO(brassdegs[keep,])
  raphdegs.RRHO = DEG2RRHO(raphdegs[keep,])
  
  return(list(brass.preRRHO = brassdegs.RRHO,
              raph.preRRHO = raphdegs.RRHO))
}

#run RRHO for cultivated vs wild progenitor
preRRHO.cultivated = prep4RRHO(degs.brass.cultivated, degs.raph.cultivated)
RRHO.brass.raph.cultivated = RRHO(preRRHO.cultivated$brass.preRRHO,
                                  preRRHO.cultivated$raph.preRRHO,
                                  stepsize = 150,
                                  labels = c("brass","raph"),
                                  alternative = "two.sided",
                                  #alternative = "enrichment",
                                  BY = T,
                                  plots = F,
                                  log10.ind = T,
                                  #outputdir = "Analysis/RNAseq/Images"
                                  )
#reverse order of heatmap columns so that they display correctly with origin at 0,0 in pheatmap
display = as.matrix(RRHO.brass.raph.cultivated$hypermat)[,order(ncol(RRHO.brass.raph.cultivated$hypermat):1)]
RRHOplot.cultivated = pheatmap(display, cluster_rows = F, cluster_cols = F, border_color = NA)
#save plot
ggsave(RRHOplot.cultivated, 
       filename = "cultivated_vs_progenitor_RRHOplot.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  42, height = 40, units = "cm")
# p.RRHO.brass.raph.cultivated = pvalRRHO(RRHO.brass.raph.cultivated, 500)
#p.RRHO.brass.raph.cultivated$pval
# save(p.RRHO.brass.raph.cultivated, file = "/home/benjamin/Documents/Brassicas_repo/Analysis/RNAseq/cultivated_vs_progenitor_RRHOpvalue.R")
# load("/home/benjamin/Documents/Brassicas_repo/Analysis/RNAseq/cultivated_vs_progenitor_RRHOpvalue.R")
# p.RRHO.brass.raph.cultivated$pval

# the p-value we get here is non-significant, but looking at the plot it seems like we might have an overlap in the downregulated genes
# let's test this by running a simple hypergeometric overlap test using just the upregulated genes and likewise for just downregulated genes
gene_overlap = function(list1, list2, background){
  
  n_A = length(list1)
  n_B = length(list2)
  n_C = length(background)
  n_A_B = length(intersect(list1,list2))
  hyp = phyper(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = FALSE)
  jac = n_A_B/(n_A+n_B-n_A_B)
  
  res = list(hypergeom = hyp,
             jaccard = jac)
  
  return(res)
  
}

#get degs and transform for blast chart
degs.brass.cultivated.up = row.names(subset(degs.brass.cultivated, padj<0.05 & log2FoldChange>0))
degs.brass.cultivated.up.sharedids = brass_raph_RBH$label[which(brass_raph_RBH$brass %in% degs.brass.cultivated.up)]
#get degs down
degs.brass.cultivated.down = row.names(subset(degs.brass.cultivated, padj<0.05 & log2FoldChange<0))
degs.brass.cultivated.down.sharedids = brass_raph_RBH$label[which(brass_raph_RBH$brass %in% degs.brass.cultivated.down)]
#get degs up
degs.raph.cultivated.up = row.names(subset(degs.raph.cultivated, padj<0.05 & log2FoldChange>0))
degs.raph.cultivated.up.sharedids = brass_raph_RBH$label[which(brass_raph_RBH$raph %in% degs.raph.cultivated.up)]
#get degs down
degs.raph.cultivated.down = row.names(subset(degs.raph.cultivated, padj<0.05 & log2FoldChange<0))
degs.raph.cultivated.down.sharedids = brass_raph_RBH$label[which(brass_raph_RBH$raph %in% degs.raph.cultivated.down)]
# no overlap in upregulated IDs
gene_overlap(degs.brass.cultivated.up.sharedids, degs.raph.cultivated.up.sharedids, brass_raph_RBH$label)
# no overlap in downregulated IDs
gene_overlap(degs.brass.cultivated.down.sharedids, degs.raph.cultivated.down.sharedids, brass_raph_RBH$label)
#so we see no overlap, but there's some weirdness here: there are far fewer Raphanus than Brass DEGs represented in the BLAST table
#could this indicate that Raph DEGs are particularly derived?
#let's first confirm that Raph DEGs are underrepresented in the table:
#get list of Raph genes that are in the table
degs.raph.cultivated.all = row.names(subset(degs.raph.cultivated,padj<0.05))
notshared.raph = row.names(degs.raph.cultivated)[which(!(row.names(degs.raph.cultivated) %in% brass_raph_RBH$raph))]
gene_overlap(notshared.raph, degs.raph.cultivated.all, row.names(degs.raph.cultivated))
#okay so raphanus DEGs are definitely overrepresented in the list of genes *without* blast best hits
#is the same for brassica?
degs.brass.cultivated.all = row.names(subset(degs.brass.cultivated,padj<0.05))
notshared.brass = row.names(degs.brass.cultivated)[which(!(row.names(degs.brass.cultivated) %in% brass_raph_RBH$brass))]
gene_overlap(notshared.brass, degs.brass.cultivated.all, row.names(degs.brass.cultivated))
#the same is not true of brassica. So perhaps we can argue that the Raphanus genes that have been selected are now more derived?
#an alternative explanation is that the relationship between R raphanistrum and R sativus is more distant than wild vs cultivated B rapa?
chiframe = data.frame(shared = !(row.names(raph.gene.counts.clean)%in%notshared.raph), 
                      DE = row.names(raph.gene.counts.clean)%in%degs.raph.cultivated.all)

plot(chisq.test(chiframe))




#run RRHO for cultivated vs wild progenitor
preRRHO.wilds = prep4RRHO(degs.brass.wilds, degs.raph.wilds)
RRHO.brass.raph.wilds = RRHO(preRRHO.wilds$brass.preRRHO,
                                  preRRHO.wilds$raph.preRRHO,
                                  stepsize = 150,
                                  labels = c("brass","raph"),
                                  alternative = "two.sided",
                                  #alternative = "enrichment",
                                  BY = T,
                                  plots = F,
                                  log10.ind = T)
#reverse order of heatmap columns so that they display correctly with origin at 0,0 in pheatmap
display = as.matrix(RRHO.brass.raph.wilds$hypermat)[,order(ncol(RRHO.brass.raph.wilds$hypermat):1)]
RRHOplot.wilds = pheatmap(display, cluster_rows = F, cluster_cols = F, border_color = NA)
#save plot
ggsave(RRHOplot.wilds, 
       filename = "progenitor_vs_wilds_RRHOplot.png",
       device = "png", path = "Analysis/RNAseq/Images/",
       width =  42, height = 40, units = "cm")
p.RRHO.brass.raph.wilds = pvalRRHO(RRHO.brass.raph.wilds, 200)
p.RRHO.brass.raph.wilds$pval
save(p.RRHO.brass.raph.wilds, file = "/home/benjamin/Documents/Brassicas_repo/Analysis/RNAseq/progenitor_vs_wilds_RRHOpvalue.R")


# out = left_join(brass_cultivated_RRHO, raph_cultivated_RRHO, "gene") %>% 
#   "colnames<-"(c("Gene","Metric1","Metric2")) %>%
#   mutate(Gene2 = Gene, Rank1 = order(Metric1), Rank2= order(Metric2))
# 
# write_delim(x=out[,c(1,4,5:6,2:3)],delim="\t", file="/home/benjamin/Downloads/test.tsv")
