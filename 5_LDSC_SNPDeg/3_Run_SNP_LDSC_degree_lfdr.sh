#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 3-00:00
#SBATCH --array=0-28
#SBATCH -p serial_requeue,xlin-lab
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mem=52Gb

module load Anaconda/5.0.1-fasrc02
source activate ldsc

declare -a LIST=("Artery_Aorta" "Artery_Coronary" "Artery_Tibial"
               "Heart_Atrial_Appendage" "Heart_Left_Ventricle" "Whole_Blood");

declare -a SUMLIST=("UKB_460K.blood_EOSINOPHIL_COUNT.sumstats.gz" "UKB_460K.blood_PLATELET_COUNT.sumstats.gz"
               "UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz" "UKB_460K.blood_RED_COUNT.sumstats.gz" "UKB_460K.blood_WHITE_COUNT.sumstats.gz"
               "PASS_HDL.sumstats.gz" "PASS_LDL.sumstats.gz");


for TISSUE in ${LIST[${SLURM_ARRAY_TASK_ID}]}
do
    for SUMSTAT in "${SUMLIST[@]}"
    do
     python /ldsc_2020/ldsc.py \
      --h2 /sumstats/$SUMSTAT \
      --ref-ld-chr /GTEX_V8/Data/LDSC_SNPDegree_Bin/$TISSUE'/bh_snp_baselineLD.',/GTEX_V8/Data/LDSC_SNPDegree_Bin/$TISSUE'/bh_snp_degree_annot.' \
      --frqfile-chr /GTEX_V8/Data/LDSC_SNPDegree_Bin/$TISSUE'/bh_snp_1000G.EUR.QC.' \
      --w-ld-chr /ldsc_2020/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
      --overlap-annot \
      --print-coefficients \
      --print-delete-vals \
      --out /GTEX_V8/Data/LDSC_SNPDegree_Bin/$TISSUE'/bh_snp_baseline.degree_annot.'$SUMSTAT
     done
done
