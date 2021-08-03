library(data.table)

#Get all tissues
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)

#Calculate degree from matrix eqtl output
calc_deg_mat <- function(eqtl_out, locat){ 
  snp_degree <- eqtl_out[,.(ind_deg=sum(indic), wgt_deg=sum(stat_abs)),by=snps]
  gene_degree <- eqtl_out[,.(ind_deg=sum(indic), wgt_deg=sum(stat_abs)),by=gene]
  names(snp_degree)[1] <- 'node'; names(gene_degree)[1] <- 'node'
  names(snp_degree)[2:3] <- paste0(names(snp_degree)[2:3],'_',locat); names(gene_degree)[2:3] <- paste0(names(gene_degree)[2:3],'_',locat)
  return(rbind(snp_degree, gene_degree))
}

#################################
##### BH-based degree
#################################
calc_deg_bh_thres <- function( eqtl_input, thres ){
  #Get the degree across all associations
  all_degree <- calc_deg_mat( eqtl_input[eqtl_input$FDR <= thres, ], 'all'); setkey(all_degree, node)
  all_cis_degree <- calc_deg_mat( eqtl_input[eqtl_input$FDR <= thres & eqtl_input$location=='cis',], 'cis'); setkey(all_cis_degree, node)
  all_trans_degree <- calc_deg_mat( eqtl_input[eqtl_input$FDR <= thres & eqtl_input$location=='trans',], 'trans'); setkey(all_trans_degree, node)
  #Now combine the degrees into a degree output 
  merge1 <- merge(all_degree, all_cis_degree, by='node', all=TRUE); degrees_all <- merge(merge1, all_trans_degree, by='node', all=TRUE)
  names(degrees_all)[2:length(names(degrees_all))] <- paste0('bh_',names(degrees_all)[2:length(names(degrees_all))])
  names(degrees_all)[2:length(names(degrees_all))] <- paste0(names(degrees_all)[2:length(names(degrees_all))], '_', thres)
  return(degrees_all)
}

#################################
##### LFDR and Qval-based degree
#################################
calc_deg_oth_thres <- function( eqtl_input, thres ){
  #Get the degree across all associations and combine the degrees into a degree output 
  #qval
  all_qval_degree <- calc_deg_mat( eqtl_input[eqtl_input$qval <= thres & (eqtl_input$qobj_cat=='cis' | eqtl_input$qobj_cat=='trans'),], 'all'); setkey(all_qval_degree, node)
  all_qval_cis_degree <- calc_deg_mat( eqtl_input[eqtl_input$qval <= thres & eqtl_input$qobj_cat=='cis',], 'cis'); setkey(all_qval_cis_degree, node)
  all_qval_trans_degree <- calc_deg_mat( eqtl_input[eqtl_input$qval <= thres & eqtl_input$qobj_cat=='trans',], 'trans'); setkey(all_qval_trans_degree, node)
  merge1 <- merge(all_qval_degree, all_qval_cis_degree, by='node', all=TRUE); degrees_qval <- merge(merge1, all_qval_trans_degree, by='node', all=TRUE)
  names(degrees_qval)[2:length(names(degrees_qval))] <- paste0('qval_',names(degrees_qval)[2:length(names(degrees_qval))])
  #lfdr
  all_lfdr_degree <- calc_deg_mat( eqtl_input[eqtl_input$lfdr <= thres & (eqtl_input$qobj_cat=='cis' | eqtl_input$qobj_cat=='trans'),], 'all'); setkey(all_lfdr_degree, node)
  all_lfdr_cis_degree <- calc_deg_mat( eqtl_input[eqtl_input$lfdr <= thres & eqtl_input$qobj_cat=='cis',], 'cis'); setkey(all_lfdr_cis_degree, node)
  all_lfdr_trans_degree <- calc_deg_mat( eqtl_input[eqtl_input$lfdr <= thres & eqtl_input$qobj_cat=='trans',], 'trans'); setkey(all_lfdr_trans_degree, node)
  merge1 <- merge(all_lfdr_degree, all_lfdr_cis_degree, by='node', all=TRUE); degrees_lfdr <- merge(merge1, all_lfdr_trans_degree, by='node', all=TRUE)
  names(degrees_lfdr)[2:length(names(degrees_lfdr))] <- paste0('lfdr_',names(degrees_lfdr)[2:length(names(degrees_lfdr))])
  #merge qval, lfdr
  degrees_all <- merge(degrees_qval, degrees_lfdr, by='node', all=TRUE)
  names(degrees_all)[2:length(names(degrees_all))] <- paste0(names(degrees_all)[2:length(names(degrees_all))], '_', thres)
  return(degrees_all)
}

#################################
##### Pnull-based degree
#################################
#Calculate degree from matrix eqtl output
calc_deg_pnull_mat <- function(eqtl_out){ 
  snp_degree <- eqtl_out[,c('snps','qobj_cat','pnull_deg')]
  deg_out <- unique(snp_degree)
  degree_out <- dcast(deg_out, snps ~ qobj_cat, value.var='pnull_deg')
  names(degree_out) <- c('node','pnull_deg_all','pnull_deg_cis','pnull_deg_trans')
  return(degree_out)
}


for (ii in 1:length(tissue_list[,1])){
  #Choose particular tissue
  tissue_id <- tissue_list[,1][ii]
  #Read in that tissues degrees for matrixeqtl, others
  eqtl_matrixeqtl <- readRDS(paste0(base.dir,"Out/Edges_matrixeqtl/Edges_matrixeqtl_",tissue_id,".Rds"))
  eqtl_matrixeqtl$stat_abs <- abs(eqtl_matrixeqtl$statistic); eqtl_matrixeqtl$indic <- 1
  char_cols <- colnames(eqtl_matrixeqtl)[which(as.vector(eqtl_matrixeqtl[,lapply(.SD, class)]) == "factor")] #make sure we don't have factors
  eqtl_matrixeqtl[,(char_cols):= lapply(.SD, as.character), .SDcols = char_cols]
  eqtl_others <- readRDS(paste0(base.dir,"Out/Edges_others/Edges_others_",tissue_id,".Rds"))
  eqtl_others$stat_abs <- abs(eqtl_others$stats); eqtl_others$indic <- 1
  #Get bh deg
  bh_deg <- Reduce(function(...) merge(..., all=TRUE),list(calc_deg_bh_thres(eqtl_matrixeqtl, 0.05), calc_deg_bh_thres(eqtl_matrixeqtl, 0.1), calc_deg_bh_thres(eqtl_matrixeqtl, 0.15),
                              calc_deg_bh_thres(eqtl_matrixeqtl, 0.2), calc_deg_bh_thres(eqtl_matrixeqtl, 0.25)))
  saveRDS(bh_deg[substr(as.character(bh_deg$node),1,3)=='chr',], file=paste0(base.dir,"Out/Degree_snp/Degree_BH_",tissue_id,".Rds"))
  saveRDS(bh_deg[substr(as.character(bh_deg$node),1,4)=='ENSG'], file=paste0(base.dir,"Out/Degree_gene/Degree_BH_",tissue_id,".Rds"))
  #Get others deg
  oth_deg <- Reduce(function(...) merge(..., all=TRUE),list(calc_deg_oth_thres(eqtl_others, 0.05), calc_deg_oth_thres(eqtl_others, 0.1), calc_deg_oth_thres(eqtl_others, 0.15),
                              calc_deg_oth_thres(eqtl_others, 0.2), calc_deg_oth_thres(eqtl_others, 0.25)))
  saveRDS(oth_deg[substr(as.character(oth_deg$node),1,3)=='chr',], file=paste0(base.dir,"Out/Degree_snp/Degree_Oth_",tissue_id,".Rds"))
  saveRDS(oth_deg[substr(as.character(oth_deg$node),1,4)=='ENSG'], file=paste0(base.dir,"Out/Degree_gene/Degree_Oth_",tissue_id,".Rds"))
  #Get pnull deg
  pnull_deg <- calc_deg_pnull_mat(eqtl_others)
  saveRDS(pnull_deg, file=paste0(base.dir,"Out/Degree_snp/Degree_Pnull_",tissue_id,".Rds"))
  saveRDS(pnull_deg, file=paste0(base.dir,"Out/Degree_gene/Degree_Pnull_",tissue_id,".Rds"))
  rm(eqtl_matrixeqtl); rm(eqtl_others); rm(bh_deg); rm(oth_deg); rm(pnull_deg); gc()
}