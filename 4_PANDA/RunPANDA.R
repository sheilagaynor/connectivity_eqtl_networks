#Take in arguments from command
indNum <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#First load in all necessary libraries
library(pandaR); library(stringr); library(data.table)

#Read in expression data
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)
tissues <- tissue_list[,1]
phenotype_data  <- fread(paste0(base.dir,"Data/GTEx_Analysis_v8_eQTL_expression_matrices/",tissues[indNum],".v8.normalized_expression.bed.gz"),data.table = F)
row.names(phenotype_data) <- gsub("\\..*","",phenotype_data$gene_id)
expression_gtex <- as.data.frame(phenotype_data[,5:ncol(phenotype_data)])

#Read in motif/ppi data
motif_gtex <- read.delim("/GTEX_V8/Data/PANDA/motif.txt", header=FALSE, stringsAsFactors=FALSE)
ppi_gtex <- read.delim("/GTEX_V8/Data/PANDA/ppi.txt", header=FALSE, stringsAsFactors=FALSE)

#Run panda
panda_gtex <- panda(motif = motif_gtex, expr = expression_gtex, ppi = ppi_gtex, 
                 remove.missing.ppi = T, remove.missing.motif = T, remove.missing.genes = T)
regTrans <- apply(panda_gtex@regNet, 2, function(x) log(exp(x) +1))
regDegree <- colSums(regTrans)

write.table(regDegree, file=paste0(base.dir,"Out/PANDA/", tissues[indNum], "_regDegree.txt"), sep="\t", col.names=FALSE, quote = FALSE)
write.table(panda_gtex@regNet, file=paste0(base.dir,"Out/PANDA/", tissues[indNum], "_regNetwork.txt"), sep="\t", quote = FALSE)
write.table(panda_gtex@coregNet, file=paste0(base.dir,"Out/PANDA/", tissues[indNum], "_coregNetwork.txt"), sep="\t", quote = FALSE)



