array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)
tissue_name <- tissue_list[,1][array_id]


for( kk in 1:5 ){
  tissue_edges <- c()
  for (ii in 1:54){
    edge_ii <- readRDS(paste0("/GTEX_V8/Out/Edges_split_others_complete/",tissue_name,"/Edges_other_",
                              tissue_name,'_test_',kk,'_',ii,'.Rds'))
    tissue_edges <- rbind(tissue_edges, edge_ii)
  }
  saveRDS(tissue_edges, paste0("/GTEX_V8/Out/Edges_split_others/",tissue_name,"/Edges_other_",
                               tissue_name,"_test_",kk,".Rds"))
  rm(tissue_edges); gc()
}


for( kk in 1:5 ){
  tissue_edges <- c()
  for (ii in 1:54){
    edge_ii <- readRDS(paste0("/GTEX_V8/Out/Edges_split_others_complete/",tissue_name,"/Edges_other_",
                              tissue_name,'_train_',kk,'_',ii,'.Rds'))
    tissue_edges <- rbind(tissue_edges, edge_ii)
  }
  saveRDS(tissue_edges, paste0("/GTEX_V8/Out/Edges_split_others/",tissue_name,"/Edges_other_",
                               tissue_name,"_train_",kk,".Rds"))
  rm(tissue_edges); gc()
}