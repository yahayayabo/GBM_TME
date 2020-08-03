#!/bin/bash -l
#SBATCH -J Epi-analysis
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=dimitrios.kyriakis@uni.lu
#SBATCH -N 1
#SBATCH --ntasks-per-node=3
#SBATCH -c 3
#SBATCH --time=0-23:59:55
#SBATCH -p batch
#SBATCH --qos=qos-batch

conda activate bioinfo_tutorial
module load swenv/default-env/devel 
module load lang/R/3.6.0-foss-2019a-bare


cd /home/users/dkyriakis/PhD/Projects/epi-scRNA/Snakemake_Analysis/


snakemake --dag | dot -Tpdf > dag.pdf
snakemake --cores 5
 
