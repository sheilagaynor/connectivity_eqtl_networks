base.dir <- "/n/holystore01/LABS/xlin/Lab/sheilagaynor/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)
tissues <- tissue_list[,1]
 
for( tissue_name in tissues ){
  tissue_edges <- c()
  for (ii in 1:54){
    edge_ii <- readRDS(paste0("/n/holystore01/LABS/xlin/Lab/sheilagaynor/GTEX_V8/Out/Edges_others_complete/Edges_other_",tissue_name,'_',ii,'.Rds'))
    tissue_edges <- rbind(tissue_edges, edge_ii)
  }
  saveRDS(tissue_edges, paste0("/n/holystore01/LABS/xlin/Lab/sheilagaynor/GTEX_V8/Out/Edges_others/Edges_others_",tissue_name,".Rds"))
  rm(tissue_edges); gc()
}
