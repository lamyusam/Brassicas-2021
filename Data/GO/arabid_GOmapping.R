#### get transcript tp GO term mapping for arabidopsis

#get data from GTF and TAIR GO file
arabid.transcript2gene =  read.delim("/home/benjamin/Documents/Brassicas_repo/Data/GO/arabid_transcript2gene.tsv", sep = " ",
                                     row.names = NULL, col.names = c("transcript","gene"))
arabid.gene2GO = read.delim("/home/benjamin/Documents/Brassicas_repo/Data/GO/arabid_gene2GO.tsv", 
                            row.names = NULL, col.names = c("gene","term"))

#for each line in transcript2gene, get GO terms in gene2GO that match that gene
arabid.GOmapping = lapply(arabid.transcript2gene$gene, function(x) arabid.gene2GO$term[which(arabid.gene2GO$gene == x)]) 
names(arabid.GOmapping) = arabid.transcript2gene$transcript

#save
save(arabid.GOmapping,file = "/home/benjamin/Documents/Brassicas_repo/Data/GO/arabid_GOmapping.Rdata")
