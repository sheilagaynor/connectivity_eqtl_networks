ND_all <- read.csv("/GTEX_V8/Out/Gene_Annot/ND_all.txt", sep="", stringsAsFactors=FALSE)
names(ND_all)[1] <- 'file_name_nd'
Taj_all <- read.csv("/GTEX_V8/Out/Gene_Annot/Taj_all.txt", sep="", stringsAsFactors=FALSE)
Taj_all$file_name_nd <- paste0('nd_', substring(Taj_all$file_name,5,nchar(Taj_all$file_name)))
gene_annot <- merge(ND_all,Taj_all,by='file_name_nd')


library(data.table); library(Matrix)

#Get tissue name
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)

corr_tab <- matrix(nrow=29,ncol=7)
for (subset_id in 1:29){
  tissue_id <- tissue_list[,1][[subset_id]]
  
  #Get the degree data for the tissue
  tissue_output_all <- readRDS(paste0(base.dir,"Out/Degree_gene/Degree_BH_",tissue_id,".Rds"))
  tissue_output_bh <- tissue_output_all[!is.na(tissue_output_all$bh_ind_deg_cis_0.05),1:5,with=F]
  rm(tissue_output_all)
  tissue_output_bh$file_name_nd <- paste0('nd_',tissue_output_bh$node,'.txt')
  gene_annot_deg <- merge(gene_annot, tissue_output_bh, by='file_name_nd')
  gene_annot_deg$transformed_deg <- log(gene_annot_deg$bh_wgt_deg_all_0.05)
  
  #Find the correlation between the annots and degrees
  corr_tab[subset_id,1] <- tissue_id
  corr_tab[subset_id,2] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg$bh_wgt_deg_all_0.05, method = "spearman")$estimate
  corr_tab[subset_id,3] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg$bh_wgt_deg_all_0.05, method = "spearman")$p.value
  corr_tab[subset_id,4] <- cor.test(gene_annot_deg$TajimaD, gene_annot_deg$bh_wgt_deg_all_0.05)$parameter + 2
  corr_tab[subset_id,5] <- cor.test(gene_annot_deg$PI, gene_annot_deg$bh_wgt_deg_all_0.05, method = "spearman")$estimate
  corr_tab[subset_id,6] <- cor.test(gene_annot_deg$PI, gene_annot_deg$bh_wgt_deg_all_0.05, method = "spearman")$p.value
  corr_tab[subset_id,7] <- cor.test(gene_annot_deg$PI, gene_annot_deg$bh_wgt_deg_all_0.05)$parameter + 2
}

saveRDS(corr_tab,'/GTEX_V8/Out/Gene_Network/Degree_NucDiv_Correlations.Rds')


library(meta)
cor_taj <- data.table(corr_tab[,c(1,2,4)], stringsAsFactors = F)
cor_taj$cor <- as.numeric(cor_taj$V2); cor_taj$n <- as.numeric(cor_taj$V3)
meta_cor_taj <- metacor(cor, n, data = cor_taj, studlab = cor_taj$V1,
                 sm = "ZCOR", method.tau = "SJ")
meta_cor_taj$TE.random
meta_cor_taj$pval.random

cor_nd <- data.table(corr_tab[,c(1,5,7)], stringsAsFactors = F)
cor_nd$cor <- as.numeric(cor_nd$V2); cor_nd$n <- as.numeric(cor_nd$V3)
meta_cor_nd <- metacor(cor, n, data = cor_nd, studlab = cor_nd$V1,
                        sm = "ZCOR", method.tau = "SJ")
meta_cor_nd$TE.random
meta_cor_nd$pval.random
