#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH --array=0-28
#SBATCH -p serial_requeue,xlin-lab,canstat-p01
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mem=32Gb

module load Anaconda/5.0.1-fasrc02
source activate ldsc

declare -a LIST=("Artery_Aorta" "Artery_Coronary" "Artery_Tibial"
               "Heart_Atrial_Appendage" "Heart_Left_Ventricle" "Whole_Blood");

for TISSUE in ${LIST[${SLURM_ARRAY_TASK_ID}]}
do
  for CHR in {1..22}
  do
    python ldsc.py \
    --l2 \
    --bfile /GTEX_V8/Data/LDSC_GeneDegree_Bin/$TISSUE'/1000G.EUR.QC.'$CHR \
    --ld-wind-cm 1 \
    --annot /GTEX_V8/Data/LDSC_GeneDegree_Bin/$TISSUE'/degree_annot.'$CHR'.annot.gz' \
    --out /GTEX_V8/Data/LDSC_GeneDegree_Bin/$TISSUE'/degree_annot.'$CHR
  done
done
