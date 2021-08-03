#Get all tissue names
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)

#Load in necessary packages
options(bitmapType='cairo')
library(WGCNA); library(data.table)

for (tissue_id in tissue_list[,1]){
#Load in the gene exp file, format correctly
phenotype_data  <- fread(paste0(base.dir,"Data/GTEx_Analysis_v8_eQTL_expression_matrices/",tissue_id,".v8.normalized_expression.bed.gz"),data.table = F)
datExpr <- as.data.frame(t(phenotype_data[,5:ncol(phenotype_data)]))

#Below taken from Peter Langfelder and Steve Horvath's
#tutorials for WGCNA package
# Choose a set of soft-thresholding powers
powers = c(c(1:14), seq(from = 16, to=34, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
png(paste0('/GTEX_V8/Out/WGCNA/Plots/',tissue_id,'_modelfit.png'))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()
}