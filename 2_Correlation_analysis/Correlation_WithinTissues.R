#Load in relevant libraries, write functions
library(data.table); library(stringr)
arrayid <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

fill_zeros = function(DT) {
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0) }
pair_corr = function(two_vec) {
  cor1 <- cor(two_vec[rowSums(two_vec)!=0,1], two_vec[rowSums(two_vec)!=0,2], use="pairwise.complete.obs", method="spearman")
  both_vecs <- which((two_vec[,1]>0)&(two_vec[,2]>0))
  cor2 <- cor(two_vec[both_vecs,1], two_vec[both_vecs,2], use="pairwise.complete.obs", method="spearman") 
  return(c(cor1,cor2))}
calculate_corr = function(test_out, train_out){
  num_degs <- ncol(test_out)-1
  joint_out <- merge(test_out, train_out, by='node', all=T)
  fill_zeros(joint_out)
  corr_out <- do.call( rbind, lapply(1:num_degs, function(x) pair_corr(joint_out[,c(x+1,x+1+num_degs),with=F]) ) )
  rownames(corr_out) <- substr(names(test_out)[2:ncol(test_out)],1,nchar(names(test_out)[2:ncol(test_out)])-5)
  return(corr_out)
}

#Get tissue names, combinations
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)
tissues <- tissue_list[,1]
tissue_name <- tissues[arrayid]

gene_deg_corr <- c(); snp_deg_corr <- c()
for (chunk in 1:5) {
  #Read in the first tissue's degree output - train
  tissue1_bh_snp <- readRDS(paste0(base.dir,"Out/Degree_split_snp/Degree_BH_",tissue_name,'_train_',chunk,".Rds"))
  setnames(tissue1_bh_snp, paste0(names(tissue1_bh_snp),'_test')); names(tissue1_bh_snp)[1] <- 'node'
  tissue1_bh_gene <- readRDS(paste0(base.dir,"Out/Degree_split_gene/Degree_BH_",tissue_name,'_train_',chunk,".Rds"))
  setnames(tissue1_bh_gene, paste0(names(tissue1_bh_gene),'_test')); names(tissue1_bh_gene)[1] <- 'node'
  #Get others deg
  tissue1_oth_snp <- readRDS(paste0(base.dir,"Out/Degree_split_snp/Degree_Oth_",tissue_name,'_train_',chunk,".Rds"))
  setnames(tissue1_oth_snp, paste0(names(tissue1_oth_snp),'_test')); names(tissue1_oth_snp)[1] <- 'node'
  tissue1_oth_gene <- readRDS(paste0(base.dir,"Out/Degree_split_gene/Degree_Oth_",tissue_name,'_train_',chunk,".Rds"))
  setnames(tissue1_oth_gene, paste0(names(tissue1_oth_gene),'_test')); names(tissue1_oth_gene)[1] <- 'node'
  #Get pnull deg
  tissue1_pnull_snp <- readRDS(paste0(base.dir,"Out/Degree_split_snp/Degree_Pnull_",tissue_name,'_train_',chunk,".Rds"))
  setnames(tissue1_pnull_snp, paste0(names(tissue1_pnull_snp),'_test')); names(tissue1_pnull_snp)[1] <- 'node'
  tissue1_pnull_gene <- readRDS(file=paste0(base.dir,"Out/Degree_split_gene/Degree_Pnull_",tissue_name,'_train_',chunk,".Rds"))
  setnames(tissue1_pnull_gene, paste0(names(tissue1_pnull_gene),'_test')); names(tissue1_pnull_gene)[1] <- 'node'
  
  #Read in the first tissue's degree output - test
  tissue2_bh_snp <- readRDS(paste0(base.dir,"Out/Degree_split_snp/Degree_BH_",tissue_name,'_test_',chunk,".Rds"))
  setnames(tissue2_bh_snp, paste0(names(tissue2_bh_snp),'_train')); names(tissue2_bh_snp)[1] <- 'node'
  tissue2_bh_gene <- readRDS(paste0(base.dir,"Out/Degree_split_gene/Degree_BH_",tissue_name,'_test_',chunk,".Rds"))
  setnames(tissue2_bh_gene, paste0(names(tissue2_bh_gene),'_train')); names(tissue2_bh_gene)[1] <- 'node'
  #Get others deg
  tissue2_oth_snp <- readRDS(paste0(base.dir,"Out/Degree_split_snp/Degree_Oth_",tissue_name,'_test_',chunk,".Rds"))
  setnames(tissue2_oth_snp, paste0(names(tissue2_oth_snp),'_train')); names(tissue2_oth_snp)[1] <- 'node'
  tissue2_oth_gene <- readRDS(paste0(base.dir,"Out/Degree_split_gene/Degree_Oth_",tissue_name,'_test_',chunk,".Rds"))
  setnames(tissue2_oth_gene, paste0(names(tissue2_oth_gene),'_train')); names(tissue2_oth_gene)[1] <- 'node'
  #Get pnull deg
  tissue2_pnull_snp <- readRDS(paste0(base.dir,"Out/Degree_split_snp/Degree_Pnull_",tissue_name,'_test_',chunk,".Rds"))
  setnames(tissue2_pnull_snp, paste0(names(tissue2_pnull_snp),'_train')); names(tissue2_pnull_snp)[1] <- 'node'
  tissue2_pnull_gene <- readRDS(file=paste0(base.dir,"Out/Degree_split_gene/Degree_Pnull_",tissue_name,'_test_',chunk,".Rds"))
  setnames(tissue2_pnull_gene, paste0(names(tissue2_pnull_gene),'_train')); names(tissue2_pnull_gene)[1] <- 'node'
  
  #Get the correlations
  snp_deg_corr_loop <- rbind( calculate_corr(tissue1_bh_snp, tissue2_bh_snp),
                              calculate_corr(tissue1_oth_snp, tissue2_oth_snp),
                              calculate_corr(tissue1_pnull_snp, tissue2_pnull_snp) )
  gene_deg_corr_loop <- rbind( calculate_corr(tissue1_bh_gene, tissue2_bh_gene),
                               calculate_corr(tissue1_oth_gene, tissue2_oth_gene),
                               calculate_corr(tissue1_pnull_gene, tissue2_pnull_gene) )
  #Name the columns
  colnames(snp_deg_corr_loop) <- c(paste0( tissue_name,'_',chunk,'_either'),
                                   paste0( tissue_name,'_',chunk,'_both') )
  colnames(gene_deg_corr_loop) <- c(paste0( tissue_name,'_',chunk,'_either'),
                                    paste0( tissue_name,'_',chunk,'_both') )
  gene_deg_corr <- cbind(gene_deg_corr,gene_deg_corr_loop)
  snp_deg_corr <- cbind(snp_deg_corr,snp_deg_corr_loop)
}

saveRDS(snp_deg_corr, paste0('/GTEX_V8/Out/Corr_within_snp/Degree_corr_within_',tissue_name,'_snp.Rds'))
saveRDS(gene_deg_corr, paste0('/GTEX_V8/Out/Corr_within_gene/Degree_corr_within_',tissue_name,'_gene.Rds'))
