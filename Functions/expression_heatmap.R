library(pheatmap)
library(RColorBrewer)

expression.heatmap = function(countdata, labels, data.phenotype = data.phenotype, main = "Untitled", pass_on = FALSE, ID_var = "Novogene_ID"){
  
  # log-transform and normalise
  heat.counts = as.matrix(countdata)
  heat.counts.vst = varianceStabilizingTransformation(heat.counts)
  heat.counts.vst.scale = t(scale(t(heat.counts.vst)))
  # scaling produces NaN values for 0-variance vectors; since these reads obviously aren't informative anyway, let's drop them
  heat.counts.vst.scale = na.omit(heat.counts.vst.scale)
  
  #get phenotypic data for heat map column annotations
  heatdata = data.frame(ID_var = colnames(countdata))
  colnames(heatdata)[1] = ID_var
  heatdata = left_join(heatdata, data.phenotype, by = ID_var)
  row.names(heatdata) = heatdata[,1]
  
  #check that supplied labels for the heatmap actually actually exist in the data.phenotype dataframe
  if(FALSE %in% (labels %in% names(heatdata)) == TRUE){
    stop(paste0("One or more invalid labels supplied. Valid labels are: ",paste(names(heatdata), collapse = ", ")))
  }
  
  #narrow down to just the phenotypic data that we want to be annotated on the plot
  heatdata = dplyr::select(heatdata, labels)
  
  #set breaks for heatmap colour scale
  breaks =  seq(-3, 3, by = 0.5)
  
  if (pass_on == FALSE) {
    plot = pheatmap(mat = heat.counts.vst.scale,
                    annotation_col = heatdata,
                    breaks = breaks,
                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
                    border_color = "white",
                    # treeheight_row = 20,
                    # treeheight_col = 15,
                    # cellwidth = 5,
                    # cellheight = 5,
                    fontsize = 5,
                    fontsize_row = 5,
                    fontsize_col = 5,
                    clustering_method = "ward.D2",
                    cluster_rows = TRUE,
                    clustering_distance_rows = "euclidean", 
                    cutree_rows = TRUE,
                    # gaps_rows = c(length(idx.gcb), length(idx.gcb) + length(idx.unc)),
                    # cluster_cols = clustCols,
                    # cutree_cols = c(2,4,6,8,10),
                    # gaps_cols = c(2,4,6,8,10),
                    clustering_distance_cols = "euclidean", 
                    legend = TRUE,
                    show_rownames = TRUE,
                    main = main)
    return(plot)
  } else if (pass_on == TRUE) {
    return(list(mat = heat.counts.vst.scale,
                annotation_col = heatdata))
  } 
}