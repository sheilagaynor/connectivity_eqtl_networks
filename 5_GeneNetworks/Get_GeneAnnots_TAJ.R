taj <- read.table('/GTEX_V8/Output/Gene_Annot/taj_files.txt', stringsAsFactors=F)
taj_genes <- substr(taj$V1, 5, nchar(taj$V1)-4)

gene_positions <- read.csv("/GTEX_V8/Output/Gene_Annot/Positions/gene_positions.tsv", sep="", stringsAsFactors=FALSE)

`%notin%` <- Negate(`%in%`)
missing_taj <- gene_positions[gene_positions$gene_id %notin% taj_genes,]
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
gene_section <- missing_taj[array_id,]

for (ii in 1:nrow(gene_section)){
  Sys.setenv(chrom=substr(gene_section[ii,1],4,nchar(gene_section[ii,1])))
  Sys.setenv(start=gene_section[ii,2])
  Sys.setenv(end=gene_section[ii,3])
  Sys.setenv(gene=gene_section[ii,4])
  system("module load vcftools/0.1.14-fasrc01; vcftools --vcf /GTEX_V8/Data/temp_selection/all_tissues.vcf --chr $chrom --from-bp $start --to-bp $end --TajimaD 2305000 --stdout > /GTEX_V8/Output/Gene_Annot/Tajimas/taj_$gene'.txt' ")
}


