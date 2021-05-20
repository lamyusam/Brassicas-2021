setwd("/home/benjamin/Documents/Brassicas_repo/")

#get phenotypic data for later
phenodata.gen2 = read.delim("Data/Phenotypic_data/morphology_data_gen2.txt", sep = ",", header = T, row.names = NULL)

#get alignment data for brassica
alignrates.brass = read.delim("Analysis/RNAseq/QC/alignment_rates/alignmentdata_brassica.tsv")
#get alignment data for raphanus
alignrates.raph = read.delim("Analysis/RNAseq/QC/alignment_rates/alignmentdata_raphanus.tsv")

#get rnaseq metadata
rna.metadata.gen2 = read.csv("/home/benjamin/Documents/Brassicas_repo/Data/RNAseq/RNASeq_sample_info.csv")[,-1]
