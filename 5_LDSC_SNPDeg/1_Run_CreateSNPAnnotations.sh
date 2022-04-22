#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -p serial_requeue,shared,xlin-lab,xlin
#SBATCH --array=1-29
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mem=20Gb

module purge
module load gcc/8.2.0-fasrc01 openmpi/3.1.1-fasrc01
module load intel-mkl/2017.2.174-fasrc01
module load R/3.6.1-fasrc01
export R_LIBS_USER=$HOME/apps/R-3.6.1-MKL:/n/helmod/apps/centos7/MPI/gcc/8.2.0-fasrc01/openmpi/3.1.1-fasrc01/R_packages/3.6.1-fasrc01

/n/home11/sgaynor/R-3.6.1/bin/R --vanilla --args $1 < /GTEX_V8/Scripts/5_LDSC/1_CreateSNPAnnotations.R
