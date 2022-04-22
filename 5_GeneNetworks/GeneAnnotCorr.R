ND_all <- read.csv("../Out/Gene_Annot/ND_all.txt", sep="", stringsAsFactors=FALSE)
names(ND_all)[1] <- 'file_name_nd'
Taj_all <- read.csv("../Out/Gene_Annot/Taj_all.txt", sep="", stringsAsFactors=FALSE)
Taj_all$file_name_nd <- paste0('nd_', substring(Taj_all$file_name,5,nchar(Taj_all$file_name)))
gene_annot <- merge(ND_all,Taj_all,by='file_name_nd')
library(data.table); library(Matrix)
#Get tissue name
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)
library(meta)

get_all_method_cor <- function(thres,zstat){
  
#BH
corr_tab <- matrix(nrow=29,ncol=7)
for (subset_id in 1:29){
  tissue_id <- tissue_list[,1][[subset_id]]
  #Get the degree data for the tissue
  tissue_output_all <- readRDS(paste0(base.dir,"Out/Degree_gene/Degree_BH_",tissue_id,".Rds"))
  tissue_output_bh <- tissue_output_all[!is.na(tissue_output_all[[paste0('bh_',zstat,'_deg_all_',thres)]]),1:27,with=F]
  rm(tissue_output_all)
  tissue_output_bh$file_name_nd <- paste0('nd_',tissue_output_bh$node,'.txt')
  gene_annot_deg <- merge(gene_annot, tissue_output_bh, by='file_name_nd')
  #Find the correlation between the annots and degrees
  corr_tab[subset_id,1] <- tissue_id
  corr_tab[subset_id,2] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg[[paste0('bh_',zstat,'_deg_all_',thres)]], method = "spearman")$estimate
  corr_tab[subset_id,3] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg[[paste0('bh_',zstat,'_deg_all_',thres)]], method = "spearman")$p.value
  corr_tab[subset_id,4] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg[[paste0('bh_',zstat,'_deg_all_',thres)]])$parameter + 2
  corr_tab[subset_id,5] <- cor.test(gene_annot_deg$PI, gene_annot_deg[[paste0('bh_',zstat,'_deg_all_',thres)]], method = "spearman")$estimate
  corr_tab[subset_id,6] <- cor.test(gene_annot_deg$PI, gene_annot_deg[[paste0('bh_',zstat,'_deg_all_',thres)]], method = "spearman")$p.value
  corr_tab[subset_id,7] <- cor.test(gene_annot_deg$PI, gene_annot_deg[[paste0('bh_',zstat,'_deg_all_',thres)]])$parameter + 2
}
cor_taj <- data.table(corr_tab[,c(1,2,4)], stringsAsFactors = F)
cor_taj$cor <- as.numeric(cor_taj$V2); cor_taj$n <- as.numeric(cor_taj$V3)
meta_cor_taj <- metacor(cor, n, data = cor_taj, studlab = cor_taj$V1,sm = "ZCOR", method.tau = "SJ")
bh_taj <- c(meta_cor_taj$TE.random,meta_cor_taj$lower.random,meta_cor_taj$upper.random)
cor_nd <- data.table(corr_tab[,c(1,5,7)], stringsAsFactors = F)
cor_nd$cor <- as.numeric(cor_nd$V2); cor_nd$n <- as.numeric(cor_nd$V3)
meta_cor_nd <- metacor(cor, n, data = cor_nd, studlab = cor_nd$V1,sm = "ZCOR", method.tau = "SJ")
bh_nd <- c(meta_cor_nd$TE.random,meta_cor_nd$lower.random,meta_cor_nd$upper.random)

#Qval
corr_tab <- matrix(nrow=29,ncol=7)
for (subset_id in 1:29){
  tissue_id <- tissue_list[,1][[subset_id]]
  #Get the degree data for the tissue
  tissue_output_all <- readRDS(paste0(base.dir,"Out/Degree_gene/Degree_Oth_",tissue_id,".Rds"))
  tissue_output_qval <- tissue_output_all[!is.na(tissue_output_all[[paste0('qval_',zstat,'_deg_all_',thres)]]),1:61,with=F]
  rm(tissue_output_all)
  tissue_output_qval$file_name_nd <- paste0('nd_',tissue_output_qval$node,'.txt')
  gene_annot_deg <- merge(gene_annot, tissue_output_qval, by='file_name_nd')
  #Find the correlation between the annots and degrees
  corr_tab[subset_id,1] <- tissue_id
  corr_tab[subset_id,2] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg[[paste0('qval_',zstat,'_deg_all_',thres)]], method = "spearman")$estimate
  corr_tab[subset_id,3] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg[[paste0('qval_',zstat,'_deg_all_',thres)]], method = "spearman")$p.value
  corr_tab[subset_id,4] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg[[paste0('qval_',zstat,'_deg_all_',thres)]])$parameter + 2
  corr_tab[subset_id,5] <- cor.test(gene_annot_deg$PI, gene_annot_deg[[paste0('qval_',zstat,'_deg_all_',thres)]], method = "spearman")$estimate
  corr_tab[subset_id,6] <- cor.test(gene_annot_deg$PI, gene_annot_deg[[paste0('qval_',zstat,'_deg_all_',thres)]], method = "spearman")$p.value
  corr_tab[subset_id,7] <- cor.test(gene_annot_deg$PI, gene_annot_deg[[paste0('qval_',zstat,'_deg_all_',thres)]])$parameter + 2
}
cor_taj <- data.table(corr_tab[,c(1,2,4)], stringsAsFactors = F)
cor_taj$cor <- as.numeric(cor_taj$V2); cor_taj$n <- as.numeric(cor_taj$V3)
meta_cor_taj <- metacor(cor, n, data = cor_taj, studlab = cor_taj$V1,sm = "ZCOR", method.tau = "SJ")
qval_taj <- c(meta_cor_taj$TE.random,meta_cor_taj$lower.random,meta_cor_taj$upper.random)
cor_nd <- data.table(corr_tab[,c(1,5,7)], stringsAsFactors = F)
cor_nd$cor <- as.numeric(cor_nd$V2); cor_nd$n <- as.numeric(cor_nd$V3)
meta_cor_nd <- metacor(cor, n, data = cor_nd, studlab = cor_nd$V1,sm = "ZCOR", method.tau = "SJ")
qval_nd <- c(meta_cor_nd$TE.random,meta_cor_nd$lower.random,meta_cor_nd$upper.random)

#FDR
corr_tab <- matrix(nrow=29,ncol=7)
for (subset_id in 1:29){
  tissue_id <- tissue_list[,1][[subset_id]]
  #Get the degree data for the tissue
  tissue_output_all <- readRDS(paste0(base.dir,"Out/Degree_gene/Degree_Oth_",tissue_id,".Rds"))
  tissue_output_lfdr <- tissue_output_all[!is.na(tissue_output_all[[paste0('lfdr_',zstat,'_deg_all_',thres)]]),1:61,with=F]
  rm(tissue_output_all)
  tissue_output_lfdr$file_name_nd <- paste0('nd_',tissue_output_lfdr$node,'.txt')
  gene_annot_deg <- merge(gene_annot, tissue_output_lfdr, by='file_name_nd')
  #Find the correlation between the annots and degrees
  corr_tab[subset_id,1] <- tissue_id
  corr_tab[subset_id,2] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg[[paste0('lfdr_',zstat,'_deg_all_',thres)]], method = "spearman")$estimate
  corr_tab[subset_id,3] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg[[paste0('lfdr_',zstat,'_deg_all_',thres)]], method = "spearman")$p.value
  corr_tab[subset_id,4] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg[[paste0('lfdr_',zstat,'_deg_all_',thres)]])$parameter + 2
  corr_tab[subset_id,5] <- cor.test(gene_annot_deg$PI, gene_annot_deg[[paste0('lfdr_',zstat,'_deg_all_',thres)]], method = "spearman")$estimate
  corr_tab[subset_id,6] <- cor.test(gene_annot_deg$PI, gene_annot_deg[[paste0('lfdr_',zstat,'_deg_all_',thres)]], method = "spearman")$p.value
  corr_tab[subset_id,7] <- cor.test(gene_annot_deg$PI, gene_annot_deg[[paste0('lfdr_',zstat,'_deg_all_',thres)]])$parameter + 2
}
cor_taj <- data.table(corr_tab[,c(1,2,4)], stringsAsFactors = F)
cor_taj$cor <- as.numeric(cor_taj$V2); cor_taj$n <- as.numeric(cor_taj$V3)
meta_cor_taj <- metacor(cor, n, data = cor_taj, studlab = cor_taj$V1,sm = "ZCOR", method.tau = "SJ")
lfdr_taj <- c(meta_cor_taj$TE.random,meta_cor_taj$lower.random,meta_cor_taj$upper.random)
cor_nd <- data.table(corr_tab[,c(1,5,7)], stringsAsFactors = F)
cor_nd$cor <- as.numeric(cor_nd$V2); cor_nd$n <- as.numeric(cor_nd$V3)
meta_cor_nd <- metacor(cor, n, data = cor_nd, studlab = cor_nd$V1,sm = "ZCOR", method.tau = "SJ")
lfdr_nd <- c(meta_cor_nd$TE.random,meta_cor_nd$lower.random,meta_cor_nd$upper.random)

return(rbind(bh_taj,bh_nd,
             qval_taj,qval_nd,
             lfdr_taj,lfdr_nd)) }

cor_05 <- get_all_method_cor('0.05','wgt')
cor_10 <- get_all_method_cor('0.1','wgt')
cor_15 <- get_all_method_cor('0.15','wgt')
cor_20 <- get_all_method_cor('0.2','wgt')
cor_25 <- get_all_method_cor('0.25','wgt')
saveRDS(rbind(cor_05,cor_10,cor_15,cor_20,cor_25),
        '../Out/Revision_Degree_GeneAnnot_Correlations_AllMethod_Wgt.Rds')

cor_05 <- get_all_method_cor('0.05','ind')
cor_10 <- get_all_method_cor('0.1','ind')
cor_15 <- get_all_method_cor('0.15','ind')
cor_20 <- get_all_method_cor('0.2','ind')
cor_25 <- get_all_method_cor('0.25','ind')
saveRDS(rbind(cor_05,cor_10,cor_15,cor_20,cor_25),
        '../Out/Revision_Degree_GeneAnnot_Correlations_AllMethod_Ind.Rds')
