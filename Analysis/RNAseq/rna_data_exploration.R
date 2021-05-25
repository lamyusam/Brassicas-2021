library(tximport)
library(DESeq2)

brass.gene.count = read.delim("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/brass.salmon.merged.gene_counts.tsv", sep = "\t", row.names = 1) 
colnames(brass.gene.count) = as.numeric(stringr::str_match(colnames(brass.gene.count),"[0-9]+"))

nextflowsheet.brass = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/nextflow_samplesheet_brass.csv") %>%
  mutate(sampleID = stringr::str_match(fastq_1,"A[0-9]+")[,1])

colnames(brass.gene.count) = nextflowsheet.brass$sampleID[match(colnames(brass.gene.count),nextflowsheet.brass$replicate)]

rna.metadata.gen2 = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info.csv")[,-1] %>%
  mutate(SpeciesShort = Species) %>%
  mutate(Domesticated = ifelse(Wild..Domesticated=="Wild","Wild","Cultivated")) %>%
  mutate(Environment = ifelse(Environment=="wheat competition","wheat","control")) %>%
  subset(substring(Species,1,3)=="Bra")

dds.gene.model = DESeqDataSetFromMatrix(countData = floor(brass.gene.count),
                                        colData = rna.metadata.gen2,
                                        design = as.formula(~Environment+Domesticated+Environment*Domesticated))

dds.gene.model.deg = DESeq(dds.gene.model, quiet = TRUE)

comparisons = resultsNames(dds.gene.model.deg)

resultstable = data.frame()

details = c()

for(i in 1:length(comparisons)){
  
  comparison = results(dds.gene.model.deg, 
                       name = comparisons[i],     
                       alpha = 0.05,
                       lfcThreshold = 0)
  
  DEGs = (length(na.omit(which(comparison$padj<0.05))))
  DEGs_pct = (length(na.omit(which(comparison$padj<0.05)))/length(na.omit(comparison$padj)))*100
  
  resultstable = rbind(resultstable, data.frame("variable" = comparisons[i],
                                                "features" = DEGs,
                                                "features_pct" = DEGs_pct))
  
  details = c(details,comparison)
  
}

resultstable

write.csv(resultstable, file = "/home/benjamin/Documents/Brassicas_repo/Analysis/RNAseq/Tables/brass_provisional_RNAseq_output.csv")

#again but with lfc threshold

resultstable = data.frame()

details = c()

for(i in 1:length(comparisons)){
  
  comparison = results(dds.gene.model.deg, 
                       name = comparisons[i],     
                       alpha = 0.05,
                       lfcThreshold = 1)
  
  DEGs = (length(na.omit(which(comparison$padj<0.05))))
  DEGs_pct = (length(na.omit(which(comparison$padj<0.05)))/length(na.omit(comparison$padj)))*100
  
  resultstable = rbind(resultstable, data.frame("variable" = comparisons[i],
                                                "features" = DEGs,
                                                "features_pct" = DEGs_pct))
  
  details = c(details,comparison)
  
}

resultstable

write.csv(resultstable, file = "/home/benjamin/Documents/Brassicas_repo/Analysis/RNAseq/Tables/brass_provisional_RNAseq_output_1lfc.csv")

