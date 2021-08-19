library(RRHO2)

raph_cultivated_RRHO2 = raph_cultivated_RRHO ; raph_cultivated_RRHO2$value[which(is.na(raph_cultivated_RRHO2$value))] = 0
brass_cultivated_RRHO2 = brass_cultivated_RRHO;  brass_cultivated_RRHO2$value[which(is.na(brass_cultivated_RRHO2$value))] = 0

RRHO2_brass_raph_cultivated = RRHO2_initialize(list1 = brass_cultivated_RRHO2,
                                              list2 = raph_cultivated_RRHO2,
                                              #stepsize = 300,
                                              #labels = c("brass","raph"),
                                              #alternative = "two.sided",
                                              #alternative = "enrichment",
                                              multipleTesting = "BH",
                                              #plots = F,
                                              log10.ind = T,
                                              #outputdir = "Analysis/RNAseq/Images"
                                              )

RRHO2_heatmap(RRHO2_brass_raph_cultivated)
dev.off()

#reverse order of heatmap columns so that they display correctly with origin at 0,0 in pheatmap
#display = as.matrix(RRHO2_brass_raph_cultivated$hypermat)[,order(ncol(RRHO2_brass_raph_cultivated$hypermat):1)]

pheatmap(as.matrix(RRHO2_brass_raph_cultivated$hypermat), cluster_rows = F, cluster_cols = F, border_color = NA)
