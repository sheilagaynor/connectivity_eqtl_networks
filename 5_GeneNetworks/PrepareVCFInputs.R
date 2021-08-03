library(Biobase)

#Get sample ids from genotype data
system("cut -d ' ' -f 2 /GTEX_V8/Data/GTEx_Analysis_v8_Genotypes/GTEX_Geno_snps_maf05_miss10_hwe1e6.fam > /GTEX_V8/Data/temp_selection/geno_samples.tsv")
geno_subjid <- read.delim("/GTEX_V8/Data/temp_selection/geno_samples.tsv",header=FALSE,stringsAsFactors=FALSE)[,1]

#Note that there are no duplicated SNPs
GTEx_bim <- read.delim("/GTEX_V8/Data/GTEx_Analysis_v8_Genotypes/GTEX_Geno_snps_maf05_miss10_hwe1e6.bim", header=FALSE, stringsAsFactors=FALSE)
length(GTEx_bim$V2) == length(unique(GTEx_bim$V2))

#Set family ID to 0 for all samples; duplicate to temp so not overwriting primary file
GTEx_fam <- read.table("/GTEX_V8/Data/temp_selection/GTEX_Geno_snps_maf05_miss10_hwe1e6.fam", quote="\"", comment.char="", stringsAsFactors=FALSE)
GTEx_fam[,1] <- 0
write.table(GTEx_fam,file="/GTEX_V8/Data/temp_selection/GTEX_Geno_snps_maf05_miss10_hwe1e6.fam",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")

#Get all of the tissue names
tissue_cov <- list.files("/GTEX_V8/Data/GTEx_Analysis_v8_eQTL_covariates")
tissues <- substr(tissue_cov, 1, nchar(tissue_cov)-18)

#Get the final snp list per tissue
tissue_keep <- c()
for(tissue_id in tissues){
  #Pull out particular tissue
  tissue_covariates <- read.delim(paste0("/GTEX_V8/Data/GTEx_Analysis_v8_eQTL_covariates/", tissue_id, ".v8.covariates.txt"), stringsAsFactors=FALSE)
  #Continue with tissues with at least 200 samples
  if( length(names(tissue_covariates))-1 > 199 ){
    tissue_keep <- c(tissue_keep, tissue_id)
    #Get rid of any duplicate samples
    remove_samples <- which(duplicated(names(tissue_covariates)))
    if(length(remove_samples) > 0){  tissue_covariates  <- tissue_covariates[,-removeSamples]  } 
    gtex_ids <- gsub("[.]", "-",  names(tissue_covariates)[2:length(names(tissue_covariates))])
    # Filter SNPs by MAF using tissue matched samples
    Sys.setenv(tissue_now=tissue_id)
    write.table(data.frame(FID=0,IID=gtex_ids),file="/GTEX_V8/Data/temp_selection/tmp_samples.tsv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")
    system("module load plink/1.90-fasrc01; plink --bfile /GTEX_V8/Data/temp_selection/GTEX_Geno_snps_maf05_miss10_hwe1e6 --const-fid --keep /GTEX_V8/Data/temp_selection/tmp_samples.tsv --geno 0.1 --maf 0.05 --hwe 1e-6 --write-snplist --out /GTEX_V8/Data/temp_selection/$tissue_now")
  } }

#Write out the final set of SNPs to keep, shared across all tissues:
snplist <- read.table(paste0("/GTEX_V8/Data/temp_selection/",tissue_keep,".snplist")[1], quote="\"", comment.char="", stringsAsFactors=FALSE)[,1]
for (tis_ind in 2:length(tissue_keep)){
  snplist <- intersect(snplist,read.table(paste0("/GTEX_V8/Data/temp_selection/",tissue_keep,".snplist")[tis_ind], quote="\"", comment.char="", stringsAsFactors=FALSE)[,1])
}
write.table(data.frame(snplist),file="/GTEX_V8/Data/temp_selection/kept_snps.snplist",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")

#Prepare final snp set
system("module load plink/1.90-fasrc01; plink --bfile /GTEX_V8/Data/temp_selection/GTEX_Geno_snps_maf05_miss10_hwe1e6 --const-fid --recode vcf --extract /GTEX_V8/Data/temp_selection/kept_snps.snplist --out /GTEX_V8/Data/temp_selection/all_tissues")


save.image("/GTEX_V8/Data/temp/QCandEqtlFiles.RData")


