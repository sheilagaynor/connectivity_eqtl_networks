#Load in necessary packages
library(WGCNA); library(data.table)
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)

#Get parameters available
tissues <- tissue_list[,1]
thresholds <- c(11,18,16,18,22,18,16,4,7,13,18,18,18,22,28,11,5,16,20,18,12,16,12,18,30,18,14,11,30)


for (indNum in 1:length(tissues)){
#Load in the gene exp file, format correctly
phenotype_data  <- fread(paste0(base.dir,"Data/GTEx_Analysis_v8_eQTL_expression_matrices/",tissues[indNum],".v8.normalized_expression.bed.gz"),data.table = F)
datExpr <- as.data.frame(t(phenotype_data[,5:ncol(phenotype_data)]))
#Calculate network
networkOut <- blockwiseModules(datExpr, power = thresholds[indNum], networkType = "signed",
                       TOMType = "signed", minModuleSize = 30, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 3)
#Calculate degree metrics
networkBlock <- intramodularConnectivity.fromExpr(datExpr, networkOut$blocks, 
                                  distFnc = "dist", distOptions = "method = 'euclidean'",
                                  networkType = "signed", power = thresholds[indNum],
                                  getWholeNetworkConnectivity = TRUE)
#Output the degree measures
write.table(cbind( phenotype_data[,4],networkBlock), file=paste0(base.dir,"/Out/WGCNA/", tissues[indNum],"_coexpDegree.txt"),
            sep="\t", row.names=FALSE, quote = FALSE)
}