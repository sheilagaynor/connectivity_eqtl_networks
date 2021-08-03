library(data.table)
#Get the list of tissue names
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)
fill_zeros = function(DT) {
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0) }
pair_corr = function(two_vec_in) {
  two_vec <- as.data.frame(two_vec_in)
  #Either
  n1 <- sum( !is.na(two_vec[,1]) & !is.na(two_vec[,2]) )
  cor1 <- cor(two_vec[rowSums(two_vec)!=0,1], two_vec[rowSums(two_vec)!=0,2], use="pairwise.complete.obs", method="spearman")
  pcor1 <- cor.test( two_vec[rowSums(two_vec)!=0,1], two_vec[rowSums(two_vec)!=0,2], method="spearman")$p.value
  #Both
  both_vecs <- which((two_vec[,1]>0)&(two_vec[,2]>0))
  n2 <- sum( !is.na(two_vec[both_vecs,1]) & !is.na(two_vec[both_vecs,2]) )
  cor2 <- cor(two_vec[both_vecs,1], two_vec[both_vecs,2], use="pairwise.complete.obs", method="spearman") 
  pcor2 <- cor.test(two_vec[both_vecs,1], two_vec[both_vecs,2], method="spearman")$p.value
  return(c(n1,cor1,pcor1,n2,cor2,pcor2))}

#For BH
degree_cor <- as.data.frame(matrix(NA, nrow=nrow(tissue_list), ncol=13))
for (ii in 1:nrow(tissue_list)){
  degree_in <- readRDS(paste0('/GTEX_V8/Out/Degree_gene/Degree_BH_',tissue_list[,1][ii],'.Rds'))
  degree_in$V1 <- gsub("\\..*","",degree_in$node)
  reg_deg <- read.delim(paste0("/GTEX_V8/Out/PANDA/",tissue_list[,1][ii],"_regDegree.txt"), 
                        sep="\t", header=FALSE, stringsAsFactors=FALSE)
  coexp_deg <- read.delim(paste0("/GTEX_V8/Out/WGCNA/",tissue_list[,1][ii],"_coexpDegree.txt"),
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)
  coexp_deg$V1 <- gsub("\\..*","",coexp_deg[,1])
  degrees_deg_panda <- merge(degree_in[,c(3,32)],reg_deg, by='V1')
  degrees_all <- merge(degrees_deg_panda,coexp_deg[,c(2,6)], by='V1')
  names(degrees_all) <- c('gene','eqtl_deg', 'panda', 'wgcna')
  setDT(degrees_all);  fill_zeros(degrees_all)
 #Get correlations
  cor_panda <- pair_corr(degrees_all[,2:3])
  cor_wgcna <- pair_corr(degrees_all[,c(2,4)])
  degree_cor[ii, 1:13] <- c(tissue_list[,1][ii], cor_panda, cor_wgcna)
}


colnames(degree_cor) <- c('tissue', 'panda_either_n', 'panda_either', 'panda_either_p', 'panda_both_n', 'panda_both', 'panda_both_p',
                          'wgcna_either_n', 'wgcna_either', 'wgcna_either_p', 'wgcna_both_n', 'wgcna_both', 'wgcna_both_p')
saveRDS(degree_cor,'/GTEX_V8/Out/Gene_Network/Degree_GeneNet_Correlations.Rds')


#Perform meta-analysis
library(meta)
panda_cor <- metacor(as.numeric(degree_cor$panda_either), as.numeric(degree_cor$panda_either_n), 
                     studlab = degree_cor$tissue, sm = "ZCOR",
                     method.tau = "SJ")
wgcna_cor <- metacor(as.numeric(degree_cor$wgcna_either), as.numeric(degree_cor$wgcna_either_n), 
                     studlab = degree_cor$tissue, sm = "ZCOR",
                     method.tau = "SJ")
