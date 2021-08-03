#Take in arrayid from command, get tissue_id
Args <- commandArgs(TRUE)
tissue_id <-as.character(Args[[1]])
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
base.dir <- "/GTEX_V8/"
set_select <- ifelse( array_id <= 269, 'test', 'train')
sample_num <- ifelse( set_select== 'test', (1 + (array_id %/% 54)), (1 + (array_id %/% 54)) - 5)
ind_num <-  (1 + (array_id %% 54)) 

#First load in all necessary libraries
library(tidyverse); library(stringr); library(magrittr)
library(data.table); library(MatrixEQTL); library(foreach) 
library(doMC); library(splines); library(MASS); library(qvalue)
registerDoMC(1)

#Make required functions for slicing and SNP positions
makeSlice_Cov <- function(fileName){
  file_in <- fread(fileName)
  covariate_sub <- file_in[-37,-1]
  sliced_data <- SlicedData$new()
  sliced_data$CreateFromMatrix(as.matrix(covariate_sub))
  rownames(sliced_data) <- file_in$ID[-37]
  return(sliced_data)}
makeSlice_Exp <- function(dataName){
  phenotype_sub <- dataName[,5:dim(dataName)[2]]
  sliced_data <- SlicedData$new()
  sliced_data$CreateFromMatrix(as.matrix(phenotype_sub))
  rownames(sliced_data) <- dataName$gene_id
  return(sliced_data)}
readSNPPos <- function(snp_pos_file){
  if(str_detect(snp_pos_file, '.gz$')){
    snp_pos_file <- str_c('zcat ', snp_pos_file)}
  snp_pos_file %>% fread %>% as.data.frame %>% set_colnames(c('snpid',	'chr', 'pos'))}

#Make function for calculating degree measures
calc_edges <- function(eqtl_out_sub, location) {
  qobj <- tryCatch(qvalue(p = eqtl_out_sub[,'pvals']), error= function(pvalues) qvalue(p =  eqtl_out_sub[,'pvals'], pi0=1))
  eqtl_out_sub[, pnull_deg := (1-qobj$pi0)] ;  eqtl_out_sub[, qval := (qobj$qvalues)];  eqtl_out_sub[, lfdr := (qobj$lfdr)]
  eqtl_out_sub[, qobj_cat := location]; eqtl_out_sub[, qobj_neqtls := nrow(eqtl_out_sub)]
  #Return the edges
  if( sum((qobj$qvalues < 0.25) | (qobj$lfdr < 0.25)) > 0 ) {
    return( eqtl_out_sub[ which((qobj$qvalues < 0.25) | (qobj$lfdr < .25)) , ] )
  } else { return( eqtl_out_sub[ which.min(qobj$qvalues$pvals) , ] ) }
}

geno_file_full  <- fread(paste0(base.dir,"Data/EQTL_Data/",tissue_id,".dosage"))
phenotype_file_full  <- fread(paste0(base.dir,"Data/GTEx_Analysis_v8_eQTL_expression_matrices/",tissue_id,".v8.normalized_expression.bed.gz"))

#Split into test and training
samp_size <- floor(0.5 * (dim(phenotype_file_full)[2]-4))
set.seed(sample_num)
obs <- sample(1:(dim(phenotype_file_full)[2]-4), size = samp_size)
#The sample set starts with 1, so may need to add up
if (set_select=="train") {  sampleSet <- (1:(dim(phenotype_file_full)[2]-4))[obs] } else {
  sampleSet <- (1:(dim(phenotype_file_full)[2]-4))[-obs] }

phenotype_data <- makeSlice_Exp(phenotype_file_full); phenotype_data$ColumnSubsample(sampleSet)
covariate_data <- makeSlice_Cov(paste0(base.dir,"Data/GTEx_Analysis_v8_eQTL_covariates/",tissue_id,".v8.covariates.txt")); covariate_data$ColumnSubsample(sampleSet)
snp_pos_data   <- readSNPPos(paste0(base.dir,"Data/EQTL_Data/",tissue_id,".pos"))
gene_pos_data <- data.frame(phenotype_file_full[,c(4,1:3)]); names(gene_pos_data)[1:2] <- c('genes','chr')
gene_pos_data$chr <- substring(gene_pos_data$chr, 4, nchar(gene_pos_data$chr))
geno_file <- geno_file_full[ seq(1,dim(snp_pos_data)[1],100000)[ind_num]:min(seq(1,dim(snp_pos_data)[1],100000)[ind_num]+99999, dim(snp_pos_data)[1]),c(1,(sampleSet+1)),with=F]
rm(geno_file_full); gc()
#Check that obs line up
table(colnames(covariate_data)==colnames(phenotype_data))
table(colnames(covariate_data)==colnames(geno_file)[2:dim(geno_file)[2]])

#Main loop to get eqtl and then prop nulls
boundIndices <- c(seq(1, dim(geno_file)[1], 1000), dim(geno_file)[1]+1) 
returnMatrix <- foreach(k=1:(length(boundIndices)-1), .combine='rbind') %dopar% {
  #Select chunk of genotpe
  geno_matrix <- as.matrix(geno_file[boundIndices[k]:(boundIndices[k+1]-1),2:(dim(geno_file)[2]),with=FALSE]) #changed bracket on boundInd
  rownames(geno_matrix) <- as.matrix(geno_file[boundIndices[k]:(boundIndices[k+1]-1),1,with=FALSE])
  genotype_data = SlicedData$new()
  genotype_data$CreateFromMatrix(geno_matrix)
  #Run matrix eQTL on slices
  matrixeqtl_obj = Matrix_eQTL_main(
    genotype_data, phenotype_data, covariate_data,
    snpspos               = snp_pos_data,      genepos               = gene_pos_data,
    output_file_name      = NULL,      output_file_name.cis  = NULL,
    pvOutputThreshold.cis = 1,      pvOutputThreshold = 1,
    min.pv.by.genesnp     = FALSE,      cisDist               = 1e6,
    useModel              = modelLINEAR,      pvalue.hist           = FALSE)
  
  eqtl_out <- data.table(snps= c(as.character(matrixeqtl_obj$trans$eqtls$snps), as.character(matrixeqtl_obj$cis$eqtls$snps)),
                         gene= c(as.character(matrixeqtl_obj$trans$eqtls$gene), as.character(matrixeqtl_obj$cis$eqtls$gene)),
                         stats= c(matrixeqtl_obj$trans$eqtls$statistic, matrixeqtl_obj$cis$eqtls$statistic),
                         pvals= c(matrixeqtl_obj$trans$eqtls$pvalue, matrixeqtl_obj$cis$eqtls$pvalue),
                         fdrs= c(matrixeqtl_obj$trans$eqtls$FDR, matrixeqtl_obj$cis$eqtls$FDR),
                         beta= c(matrixeqtl_obj$trans$eqtls$FDR, matrixeqtl_obj$cis$eqtls$beta),
                         loc= c(rep("trans", length(matrixeqtl_obj$trans$eqtls$pvalue)),rep("cis",length(matrixeqtl_obj$cis$eqtls$pvalue))))
  setkey(eqtl_out, snps);  snpList <- unique(eqtl_out[,snps])
  
  foreach(ii=1:length(snpList), .combine='rbind') %do% { 
    all_edge <- calc_edges( eqtl_out[.(snpList[ii]),] , 'all') 
    trans_edge <- calc_edges( eqtl_out[.(snpList[ii]),][eqtl_out[.(snpList[ii]),loc]=='trans'] , 'trans') 
    if (length(which(eqtl_out[.(snpList[ii]),loc]=='cis'))>1){
      cis_edge <- calc_edges( eqtl_out[.(snpList[ii]),][eqtl_out[.(snpList[ii]),loc]=='cis'] , 'cis')
      rbind( all_edge, trans_edge, cis_edge) } else {
        rbind( all_edge, trans_edge) }
  } 
}

saveRDS(returnMatrix, file=paste0(base.dir,"Out/Edges_split_others_complete/",tissue_id,"/Edges_other_",
                                  tissue_id, '_', set_select, '_', sample_num, '_', ind_num, ".Rds"))
