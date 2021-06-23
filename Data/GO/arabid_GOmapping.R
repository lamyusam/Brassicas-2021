#### get transcript tp GO term mapping for arabidopsis

setwd("/home/benjamin/Documents/Brassicas_repo")

#get data from GTF and TAIR GO file
arabid.transcript2gene =  read.delim("Data/GO/arabid_transcript2gene.tsv", sep = " ",
                                     row.names = NULL, col.names = c("transcript","gene"))
arabid.gene2GO = read.delim("Data/GO/arabid_gene2GO.tsv", 
                            row.names = NULL, col.names = c("gene","term"))

#for each line in transcript2gene, get GO terms in gene2GO that match that gene
arabid.GOmapping = lapply(arabid.transcript2gene$gene, function(x) arabid.gene2GO$term[which(arabid.gene2GO$gene == x)]) 
names(arabid.GOmapping) = arabid.transcript2gene$transcript

#save
save(arabid.GOmapping,file = "Data/GO/arabid_GOmapping.Rdata")

load("Data/GO/arabid_GOmapping.Rdata")
#import gene2LOC adapted from raphanus gff
raph.exon2LOC = read.delim("Data/RNAseq/raph_exon2LOC.tsv", col.names = c("exon","gene"))

#### now attach Raphanus data
raph.arabid.BLAST = read.delim("Data/GO/arabid_raph_RBH.txt", sep = " ", col.names = c("raphanus","arabid"))
#reduce to short transcript identifiers
raph.arabid.BLAST = mutate(raph.arabid.BLAST, 
                           arabid = str_match(raph.arabid.BLAST$arabid,"(NM_[^_\n]+)|([rt]rna_[^_\n]+)|(mrna_[1-9]+)")[,1],
                           raphanus = str_match(raph.arabid.BLAST$raphanus,"(X[MR]_[^_\n]+)|([rt]rna_[^_\n]+)|(mrna_[1-9]+)")[,1])

#function to apply over the list of gene IDs, attaching arabidopsis GO terms where possible
lapfun = function(x){
  unname(unlist(arabid.GOmapping[raph.arabid.BLAST$arabid[which(raph.arabid.BLAST$raphanus %in% 
                                                  raph.exon2LOC$exon[which(raph.exon2LOC$gene==x)])]]))
}

#apply function over list to generate GO mapping object
GOmapping.raph = lapply(rownames(raph.gene.counts), lapfun)
names(GOmapping.raph) = rownames(raph.gene.counts)


# raph.arabid.BLAST
# "NM_001198228.1" %in% names(arabid.GOmapping)
# arabid.GOmapping["NM_001198228.1"]
