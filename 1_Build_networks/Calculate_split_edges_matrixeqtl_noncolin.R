#Take in arrayid from command, get tissue_id
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)
tissue_id <- tissue_list[(1 + (array_id %/% 10)),]
set_select <- ifelse( (1 + (array_id %% 10)) <= 5, 'test', 'train')
sample_num <- ifelse( set_select=='test', (1 + (array_id %% 10)), (1 + (array_id %% 10)) - 5)

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

geno_file_full  <- fread(paste0(base.dir,"Data/EQTL_Data/",tissue_id,".dosage"))
phenotype_file_full  <- fread(paste0(base.dir,"Data/GTEx_Analysis_v8_eQTL_expression_matrices/",tissue_id,".v8.normalized_expression.bed.gz"))

#Split into test and training
samp_size <- floor(0.5 * (dim(phenotype_file_full)[2]-4))
set.seed(sample_num)
obs <- sample(1:(dim(phenotype_file_full)[2]-4), size = samp_size)
#The sample set starts with 1, so may need to add up
if (set_select=="train") {  sample_set <- (1:(dim(phenotype_file_full)[2]-4))[obs] } else {
  sample_set <- (1:(dim(phenotype_file_full)[2]-4))[-obs] }

phenotype_data <- makeSlice_Exp(phenotype_file_full); phenotype_data$ColumnSubsample(sample_set)
covariate_data <- makeSlice_Cov(paste0(base.dir,"Data/GTEx_Analysis_v8_eQTL_covariates/",tissue_id,".v8.covariates.txt")); covariate_data$ColumnSubsample(sample_set)
snp_pos_data   <- readSNPPos(paste0(base.dir,"Data/EQTL_Data/",tissue_id,".pos"))
gene_pos_data <- data.frame(phenotype_file_full[,c(4,1:3)]); names(gene_pos_data)[1:2] <- c('genes','chr')
gene_pos_data$chr <- substring(gene_pos_data$chr, 4, nchar(gene_pos_data$chr))
geno_file <- geno_file_full[ ,c(1,(sample_set+1)),with=F] 
rm(geno_file_full); gc()
#Check that obs line up
table(colnames(covariate_data)==colnames(phenotype_data))
table(colnames(covariate_data)==colnames(geno_file)[2:dim(geno_file)[2]])


  #Prepare genotpe
  geno_matrix <- as.matrix(geno_file[,2:(dim(geno_file)[2]),with=FALSE]) 
  rownames(geno_matrix) <- as.matrix(geno_file[,1,with=FALSE]) 
  rm(geno_file); gc()
  genotype_data = SlicedData$new()
  genotype_data$CreateFromMatrix(geno_matrix)
  genotype_data$ResliceCombined(sliceSize = 5000L)
  rm(geno_matrix); gc()
  #Check that obs line up
  table(colnames(covariate_data)==colnames(genotype_data))
  
  
  #Run matrix eQTL on slices
  matrixeqtl_obj = Matrix_eQTL_main(
    genotype_data, phenotype_data, covariate_data,
    snpspos               = snp_pos_data,      genepos               = gene_pos_data,
    output_file_name      = NULL,      output_file_name.cis  = NULL,
    pvOutputThreshold.cis = 1e-1,      pvOutputThreshold = 1e-4,
    min.pv.by.genesnp     = FALSE,      cisDist               = 1e6,
    useModel              = modelLINEAR,      pvalue.hist           = FALSE )
  
  trans <- cbind(matrixeqtl_obj$trans$eqtls, 'trans');  cis <- cbind(matrixeqtl_obj$cis$eqtls, 'cis')
  names(trans)[7] <- 'location'; names(cis)[7] <- 'location'
  eqtl_out <- rbind(trans, cis); rm(matrixeqtl_obj); rm(trans); rm(cis)
  setDT(eqtl_out)
  eqtls <- eqtl_out[eqtl_out$FDR < 0.25,]; rm(eqtl_out)
  
  saveRDS(eqtls, file=paste0(base.dir,"Out/Edges_split_matrixeqtl/Edges_split_matrixeqtl_",
                             tissue_id, '_', sample_num, '_', set_select, ".Rds"))
  