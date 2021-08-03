###############################################
#### Step 1: Create annot files
###############################################
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
readSNPPos <- function(snp_pos_file){
  if(str_detect(snp_pos_file, '.gz$')){
    snp_pos_file <- str_c('zcat ', snp_pos_file)}
  snp_pos_file %>% fread %>% as.data.frame %>% set_colnames(c('node',	'chr', 'pos'))}

################
#In R
library(data.table);library(tidyverse); library(stringr); library(magrittr); library(scales)
#Get the lifted positions
bed_map <- fread('/GTEX_V8/Data/GTEx_Analysis_v8_Genotypes/GTEX_Geno_snps_maf05_miss10_hwe1e6.mapping')
names(bed_map)[which(names(bed_map)=='variant_id')] <- c('node')
setkey(bed_map,node)
#Get the list of tissue names
base.dir <- "/GTEX_V8/"
tissue_list <- read.table(paste0(base.dir,"Data/EQTL_Data_All/tissue_list.txt"), quote="\"", comment.char="", stringsAsFactors=FALSE)
#Loop through to prepare .annot for each tissue
tissue_id <- tissue_list[,1][array_id]
Sys.setenv(tissue_id=tissue_id)
degree_vals <- readRDS(paste0('/GTEX_V8/Out/Degree_snp/Degree_BH_',tissue_id,'.Rds'))
degree_vals$network_annot <- as.numeric(degree_vals$bh_ind_deg_all_0.05 >= 0)
degree_vals$degree_annot <- as.numeric(degree_vals$bh_wgt_deg_all_0.05 >= quantile(degree_vals$bh_wgt_deg_all_0.05, 0.75, na.rm=T))
setkey(degree_vals,node)
degree_mapped <- merge(degree_vals[,c('node','network_annot','degree_annot')], bed_map[,c('node','chr','variant_id_b37')], by='node')

#Make annotations for each chromosome
  for (chr in 1:22){
    Sys.setenv(chr=chr)
    #Prepare annots
    chr_ind <- which(substring(degree_mapped$chr,4,nchar(degree_mapped$chr)) == chr)
    degree_chr <- degree_mapped[chr_ind,]
    bim_in <- fread(paste0('/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.',chr,'.bim'))
    names(bim_in) <- c('CHR','SNP','CM','BP','REF','ALT')
    degree_chr[, c("CHROM", "BP_STR", 'REF', 'ALT', 'BUILD') := tstrsplit(variant_id_b37, "_", fixed=TRUE)]
    degree_chr$BP <- as.numeric(degree_chr$BP_STR)
    key_col <- c('BP')
    #Remove duplicates
    bed_chr_dup <- merge(bim_in, degree_chr[,c('network_annot','degree_annot','BP')], by=key_col, all.x=T)
    set(bed_chr_dup,which(is.na(bed_chr_dup[,'network_annot'])),7L,0)
    set(bed_chr_dup,which(is.na(bed_chr_dup[,'degree_annot'])),8L,0)
    bed_chr <- bed_chr_dup[,.SD[which.max(degree_annot)],by=SNP]
    #Get relevant annot files and save
    baseld <- fread( paste0('/ldsc/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.',chr,'.annot.gz'))
    write.table(baseld[baseld$SNP %in% bed_chr$SNP,],
                file=gzfile(paste0('/GTEX_V8/Data/LDSC_SNPDegree_Bin/',tissue_id,'/baselineLD.',chr,'.annot.gz')),
                row.names = FALSE, sep="\t", quote = FALSE)
    write.table(bed_chr[,c('CHR','BP','SNP','CM','network_annot')],
                file=gzfile(paste0('/GTEX_V8/Data/LDSC_SNPDegree_Bin/',tissue_id,'/network_annot.',chr,'.annot.gz')),
                row.names = FALSE, sep="\t", quote = FALSE)
    write.table(bed_chr[,c('CHR','BP','SNP','CM','degree_annot')],
                file=gzfile(paste0('/GTEX_V8/Data/LDSC_SNPDegree_Bin/',tissue_id,'/degree_annot.',chr,'.annot.gz')),
                row.names = FALSE, sep="\t", quote = FALSE)
    #Prepare frq, bim/bed/fam
    snp_list <- fread(paste0('/GTEX_V8/Data/LDSC_SNPDegree_Bin/',tissue_id,'/network_annot.',chr,'.annot.gz'))
    write.table(data.frame(snp_list$SNP),file=paste0("/GTEX_V8/Data/LDSC_SNPDegree_Bin/",tissue_id,"/kept_snps.snplist"),row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")
    frq <- fread(paste0('/ldsc/1000G_Phase3_frq/1000G.EUR.QC.',chr,'.frq'))
    write.table(frq[frq$SNP %in% snp_list$SNP,],
                file=(paste0('/GTEX_V8/Data/LDSC_SNPDegree_Bin/',tissue_id,'/1000G.EUR.QC.',chr,'.frq')),
                row.names = FALSE, sep="\t", quote = FALSE)
    #Subset plink files by tissue/chrom   
    system("module load plink/1.90-fasrc01; plink --bfile /ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr --extract /GTEX_V8/Data/LDSC_SNPDegree_Bin/$tissue_id'/kept_snps.snplist' --make-bed --out /GTEX_V8/Data/LDSC_SNPDegree_Bin/$tissue_id'/1000G.EUR.QC.'$chr ")
  }
