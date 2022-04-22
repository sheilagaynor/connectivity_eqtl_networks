###############################################
#### Step 1: Create annot files
###############################################
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
args <- commandArgs(TRUE)
meas <- as.character(args[[1]]) 

################
#In R
library(data.table); library(scales)
#Get the lifted positions
bed_map <- fread('/GTEX_V8/Data/GTEx_Analysis_v8_eQTL_expression_matrices/lifted_gtex_genes.bed')
names(bed_map) <- c( 'CHR', 'START', 'END', 'node')
bed_map$START_WINDOW <- bed_map$START - 25000 
bed_map$END_WINDOW <- bed_map$END + 25000
setkey(bed_map,node)
#Get the list of tissue names
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)
#Loop through to prepare .annot for each tissue
tissue_id <- tissue_list[,1][array_id]
Sys.setenv(tissue_id=tissue_id)
Sys.setenv(meas=meas)
if (meas == 'bh'){
  degree_in <- readRDS(paste0('/GTEX_V8/Out/Degree_gene/Degree_BH_',tissue_id,'.Rds'))
  degree_in$network_annot <- as.numeric(degree_in$bh_ind_deg_all_0.05 >= 0)
  degree_in$degree_annot <- as.numeric(degree_in$bh_wgt_deg_all_0.05 >= quantile(degree_in$bh_wgt_deg_all_0.05, 0.75, na.rm=T))
} else if (meas == 'lfdr'){
  degree_in <- readRDS(paste0('/GTEX_V8/Out/Degree_gene/Degree_Oth_',tissue_id,'.Rds'))
  degree_in$network_annot <- as.numeric(degree_in$lfdr_ind_deg_all_0.05 >= 0)
  degree_in$degree_annot <- as.numeric(degree_in$lfdr_wgt_deg_all_0.05 >= quantile(degree_in$lfdr_wgt_deg_all_0.05, 0.75, na.rm=T))
} else if (meas == "qval"){
  degree_in <- readRDS(paste0('/GTEX_V8/Out/Degree_gene/Degree_Oth_',tissue_id,'.Rds'))
  degree_in$network_annot <- as.numeric(degree_in$qval_ind_deg_all_0.05 >= 0)
  degree_in$degree_annot <- as.numeric(degree_in$qval_wgt_deg_all_0.05 >= quantile(degree_in$qval_wgt_deg_all_0.05, 0.75, na.rm=T))
}
setkey(degree_in,node)
degree_mapped <- merge(degree_in[,c('node','network_annot','degree_annot')], bed_map, by='node', all.y=T)

  #Make annotations for each chromosome
  for (chr in 1:22){
    Sys.setenv(chr=chr)
    #Prepare annots
    chr_ind <- which(substring(degree_mapped$CHR,4,nchar(degree_mapped$CHR)) == chr)
    degree_chr <- degree_mapped[chr_ind,]
    bim_in <- fread(paste0('/ldsc_2020/1000G_EUR_Phase3_plink/1000G.EUR.QC.',chr,'.bim'))
    names(bim_in) <- c('CHR','SNP','CM','BP','REF','ALT')
    bim_in$START_WINDOW <- bim_in$BP;   bim_in$END_WINDOW <- bim_in$BP
    key_col <- c('START_WINDOW', 'END_WINDOW')
    setkeyv(degree_chr, key_col); setkeyv(bim_in, key_col)
    bed_chr_over <- foverlaps(bim_in, degree_chr, type="any", nomatch=NA)
    bed_chr_over <- bed_chr_over[!is.na(bed_chr_over$node) & !is.na(bed_chr_over$network_annot),]
    #Remove duplicates
    setkey(bed_chr_over, 'SNP')
    bed_no_dup <- bed_chr_over[,.SD[which.max(degree_annot)],by=SNP]
    setkey(bed_no_dup, 'SNP'); setkey(bim_in, 'SNP')
    bed_chr <- merge(bim_in, bed_no_dup[,c('network_annot','degree_annot','SNP')], by='SNP', all.x=T)
    set(bed_chr,which(is.na(bed_chr[,'network_annot'])),9L,0)
    set(bed_chr,which(is.na(bed_chr[,'degree_annot'])),10L,0)
    bed_chr <- setorder(bed_chr, BP)
    #Get relevant annot files and save
    baseld <- fread( paste0('/ldsc_2020/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.',chr,'.annot.gz'))
    write.table(baseld[baseld$SNP %in% bed_chr$SNP,],
                file=gzfile(paste0('/GTEX_V8/Data/LDSC_GeneDegree_Bin/',tissue_id,'/',meas,'_gene_baselineLD.',chr,'.annot.gz')),
                row.names = FALSE, sep="\t", quote = FALSE)
    write.table(bed_chr[,c('CHR','BP','SNP','CM','network_annot')],
                file=gzfile(paste0('/GTEX_V8/Data/LDSC_GeneDegree_Bin/',tissue_id,'/',meas,'_gene_network_annot.',chr,'.annot.gz')),
                row.names = FALSE, sep="\t", quote = FALSE)
    write.table(bed_chr[,c('CHR','BP','SNP','CM','degree_annot')],
                file=gzfile(paste0('/GTEX_V8/Data/LDSC_GeneDegree_Bin/',tissue_id,'/',meas,'_gene_degree_annot.',chr,'.annot.gz')),
                row.names = FALSE, sep="\t", quote = FALSE)
    #Prepare frq, bim/bed/fam
    snp_list <- fread(paste0('/GTEX_V8/Data/LDSC_GeneDegree_Bin/',tissue_id,'/',meas,'_gene_network_annot.',chr,'.annot.gz'))
    write.table(data.frame(snp_list$SNP),file=paste0("/GTEX_V8/Data/LDSC_GeneDegree_Bin/",tissue_id,'/',meas,"_gene_kept_snps.snplist"),row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")
    frq <- fread(paste0('/ldsc_2020/1000G_Phase3_frq/1000G.EUR.QC.',chr,'.frq'))
    write.table(frq[frq$SNP %in% snp_list$SNP,],
                file=(paste0('/GTEX_V8/Data/LDSC_GeneDegree_Bin/',tissue_id,'/',meas,'_gene_1000G.EUR.QC.',chr,'.frq')),
                row.names = FALSE, sep="\t", quote = FALSE)
    #Subset plink files by tissue/chrom   
    system("module load plink/1.90-fasrc01; plink --bfile /ldsc_2020/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr \\
           --extract /GTEX_V8/Data/LDSC_GeneDegree_Bin/$tissue_id'/'$meas'_gene_kept_snps.snplist' \\
           --make-bed --out /GTEX_V8/Data/LDSC_GeneDegree_Bin/$tissue_id'/'$meas'_gene_1000G.EUR.QC.'$chr ")
  }
