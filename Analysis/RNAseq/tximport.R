library("tximport")
library("GenomicFeatures")
library("tidyverse")

setwd("/home/benjamin/Documents/Brassicas_repo")

#make tximport databse for raphanus, matching transcripts to genes in gff
txdb = makeTxDbFromGFF("Data/RNAseq/Raph.gff")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
#check that it worked as expected
head(tx2gene)
length(unique(tx2gene$GENEID))
length(unique(tx2gene$TXNAME))
#for some reason tximport removes prefixes that are needed for matching to the salmon files, so add these back
#we also have to re-order the columns for tximport to read correctly (no idea why this isn't automatic)
tx2gene = mutate(na.omit(tx2gene), 
                 GENEID = paste0("gene-",GENEID),
                 TXNAME = paste0("rna-",TXNAME))[,c(2,1)]

#now we need a list of Salmon transcript output files
salmfiles.raph = list.files("Data/RNAseq/raph_quant", full.names = T)
names(salmfiles.raph) = str_match(salmfiles.raph,"R[0-9]+") #set names
txi.raph = tximport(salmfiles.raph, type = "salmon", tx2gene = tx2gene)
#if passing to DESeq2, the recommendation is to use DESeqDataSetFromTximport to round the estimated values...
#however, per Mike Love, all DESeq2 does is round the values anyway, so there's no harm in doing it the old-fashioned way
txi.raph.counts = round(txi.raph$counts)
#convert nextflow labels to original metadata labels
nextflowsheet.raph = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/nextflow_samplesheet_raph.csv") %>%
  mutate(sampleID = stringr::str_match(fastq_1,"A[0-9]+")[,1],
         replicate = paste0("R",replicate))
colnames(txi.raph.counts) = nextflowsheet.raph$sampleID[match(colnames(txi.raph.counts),nextflowsheet.raph$replicate)]

#save
write.csv(txi.raph.counts,file = "Data/RNAseq/raph.gene.counts.csv")





####### brapa

#make tximport databse for brassica, matching transcripts to genes in gff
txdb = makeTxDbFromGFF("Data/RNAseq/oldBrapa.gff")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
#check that it worked as expected
head(tx2gene)
length(unique(tx2gene$GENEID))
length(unique(tx2gene$TXNAME))
#for some reason tximport removes prefixes that are needed for matching to the salmon files, so add these back
#we also have to re-order the columns for tximport to read correctly (no idea why this isn't automatic)
tx2gene = mutate(na.omit(tx2gene), 
                 GENEID = paste0("gene-",GENEID),
                 TXNAME = paste0("rna-",TXNAME))[,c(2,1)]

#now we need a list of Salmon transcript output files
salmfiles.brapa = list.files("Data/RNAseq/brapa_quant/", full.names = T)
names(salmfiles.brapa) = str_match(salmfiles.brapa,"R[0-9]+") #set names
txi.brapa = tximport(salmfiles.brapa, type = "salmon", tx2gene = tx2gene)
#if passing to DESeq2, the recommendation is to use DESeqDataSetFromTximport to round the estimated values...
#however, per Mike Love, all DESeq2 does is round the values anyway, so there's no harm in doing it the old-fashioned way
txi.brapa.counts = round(txi.brapa$counts)
#convert nextflow labels to original metadata labels
nextflowsheet.brapa = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/nextflow_samplesheet_brass.csv") %>%
  mutate(sampleID = stringr::str_match(fastq_1,"A[0-9]+")[,1],
         replicate = paste0("R",replicate))
colnames(txi.brapa.counts) = nextflowsheet.brapa$sampleID[match(colnames(txi.brapa.counts),nextflowsheet.brapa$replicate)]

#save
#write.csv(txi.brapa.counts,file = "Data/RNAseq/brapa.gene.counts.csv")
