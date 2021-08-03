nd <- read.table('/GTEX_V8/Output/Gene_Annot/nd_files.txt', stringsAsFactors=F)
nd_genes <- substr(nd$V1, 4, nchar(nd$V1)-4)

gene_positions <- read.csv("/GTEX_V8/Output/Gene_Annot/Positions/gene_positions.tsv", sep="", stringsAsFactors=FALSE)

`%notin%` <- Negate(`%in%`)
missing_nd <-  gene_positions[gene_positions$gene_id  %notin% nd_genes,]
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
gene_section <- missing_nd[array_id,]

for (ii in 1:nrow(gene_section)){
  Sys.setenv(chrom=substr(gene_section[ii,1],4,nchar(gene_section[ii,1])))
  Sys.setenv(start=gene_section[ii,2])
  Sys.setenv(end=gene_section[ii,3])
  Sys.setenv(gene=gene_section[ii,4])
  system("module load vcftools/0.1.14-fasrc01; vcftools --vcf /GTEX_V8/Data/temp_selection/all_tissues.vcf --chr $chrom --from-bp $start --to-bp $end --window-pi 2305000 --stdout > /GTEX_V8/Output/Gene_Annot/ND/nd_$gene'.txt' ")
}


