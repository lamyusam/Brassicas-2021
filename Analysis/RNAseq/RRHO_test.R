library(RRHO)
library(stringr)

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

# degs.brass.cultivated = results(brassrapa.dds, name="treatment_Control_vs_Wheat", alpha = 0.05)     
# degs.raph.cultivated = results(raphsubset.dds, name="treatment_Control_vs_Wheat", alpha = 0.05) 

# degs.brass.cultivated = results(brasswilds.dds, name="wild.ancestorTRUE", alpha = 0.05)     
# degs.raph.cultivated = results(raphwilds.dds, name="wild.ancestorTRUE", alpha = 0.05) 

# degs.brass.cultivated = results(brasswilds.dds, name="treatment_Control_vs_Wheat", alpha = 0.05)     
# degs.raph.cultivated = results(raphwilds.dds, name="treatment_Control_vs_Wheat", alpha = 0.05) 

#remove rows that aren't in the matching table
degs.brass.cultivated = subset(degs.brass.cultivated, row.names(degs.brass.cultivated) %in% brass_raph_RBH$brass)
#relabel each output's rownames to match the labels in the Blast output table
row.names(degs.brass.cultivated) = brass_raph_RBH$label[match(brass_raph_RBH$brass, row.names(degs.brass.cultivated))]

#remove rows that aren't in the matching table
degs.raph.cultivated = subset(degs.raph.cultivated, row.names(degs.raph.cultivated) %in% brass_raph_RBH$raph)
#relabel each output's rownames to match the labels in the Blast output table
row.names(degs.raph.cultivated) = brass_raph_RBH$label[match(brass_raph_RBH$raph, row.names(degs.raph.cultivated))]

#reorder the frames to match
degs.raph.cultivated = degs.raph.cultivated[row.names(degs.brass.cultivated),]

#now we should be able to input to RRHO

#before transforming for RRHO, we might wish to remove the set of genes that have p values close to 1 in both comparisons
#this may improve the analysis because genes with very low significance are unlikely to rank meaningfully
#combine two pvalue lists
both = cbind(degs.brass.cultivated$pvalue, degs.raph.cultivated$pvalue) 
#remove rows that have p-value above threshold in both comparisons
keep =  row.names(degs.brass.cultivated)[(rowSums(both > 0.5 | is.na(both))!=2)]
print(paste0("Removing ",nrow(degs.brass.cultivated)-length(keep)," genes with low significance in both comparisons."))

#RRHO prep function
DEG2RRHO = function(contrast){
  
  RRHO = data.frame(gene = row.names(contrast),
                    value = -log10(contrast$pvalue)*sign(contrast$log2FoldChange))
                    #value = contrast$log2FoldChange)
  return(RRHO)
  
}

brass_cultivated_RRHO = DEG2RRHO(degs.brass.cultivated[keep,])
raph_cultivated_RRHO = DEG2RRHO(degs.raph.cultivated[keep,])

RRHO_brass_raph_cultivated = RRHO(brass_cultivated_RRHO,
                                  raph_cultivated_RRHO,
                                  stepsize = 200,
                                  labels = c("brass","raph"),
                                  alternative = "two.sided",
                                  #alternative = "enrichment",
                                  BY = T,
                                  plots = F,
                                  log10.ind = T,
                                  #outputdir = "Analysis/RNAseq/Images"
                                  )


#lattice::levelplot(RRHO_brass_raph_cultivated$hypermat)
# col = colorRampPalette(brewer.pal(8, "RdYlBu"))(25)
# heatmap(RRHO_brass_raph_cultivated$hypermat, Rowv = NA, Colv = NA, col = col, labRow = NA, labCol = NA)

#reverse order of heatmap columns so that they display correctly with origin at 0,0 in pheatmap
display = as.matrix(RRHO_brass_raph_cultivated$hypermat)[,order(ncol(RRHO_brass_raph_cultivated$hypermat):1)]

pheatmap(display, cluster_rows = F, cluster_cols = F, border_color = NA)

p_RRHO_brass_raph_cultivated = pvalRRHO(RRHO_brass_raph_cultivated, 1)
p_RRHO_brass_raph_cultivated$pval

out = left_join(brass_cultivated_RRHO, raph_cultivated_RRHO, "gene") %>% 
  "colnames<-"(c("Gene","Metric1","Metric2")) %>%
  mutate(Gene2 = Gene, Rank1 = order(Metric1), Rank2= order(Metric2))

write_delim(x=out[,c(1,4,5:6,2:3)],delim="\t", file="/home/benjamin/Downloads/test.tsv")
