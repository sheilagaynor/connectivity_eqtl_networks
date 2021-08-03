#Read in gencode gene information
gencode.v26.annotation <- read.delim("/Data/gencode.v26.annotation.gtf", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
gencode <- gencode.v26.annotation[gencode.v26.annotation$V3 =='gene',]
rm(gencode.v26.annotation)
gencode$gene_id_labeled <- gsub(";.*","",gencode$V9)
gencode$gene_id <- substr(gencode$gene_id_labeled, 9, nchar(gencode$gene_id_labeled))
#Write out the bed data needed for liftover
bed_format <- data.frame( gencode$V1, gencode$V4, gencode$V5, gencode$gene_id)
write.table(bed_format, "/Data/gtex_genes.bed", col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)

system("./liftOver gtex_genes.bed hg38ToHg19.over.chain lifted_gtex_genes.bed unlifted.bed")
