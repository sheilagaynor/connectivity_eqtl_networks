library(data.table)
gencode <- fread('gencode.v26.annotation.gtf')
gencode_genes <- gencode[gencode$V3=='gene',]
text_col <- strsplit(gencode_genes$V9, ';')
gene_id <- unlist(lapply(text_col, function(l) l[[1]]))
gene_id_trimmed <- substr(gene_id, 10, nchar(gene_id)-1)
gene_pos <- cbind(gencode_genes[,c(1,4:5)], gene_id_trimmed)
names(gene_pos) <- c('chr', 'start', 'end', 'gene_id')
gene_pos$width <- gene_pos$end - gene_pos$start
write.table(gene_pos, 'gene_positions.tsv', quote = F, row.names = F)

