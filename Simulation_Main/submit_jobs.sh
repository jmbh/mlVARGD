#!/bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=thin
#SBATCH --time=13:00:00

module purge
module load 2021 R/4.1.0-foss-2021a

cp -r "$HOME"/mlVARGD_Sim13 "$TMPDIR"
cd "$TMPDIR"/mlVARGD_Sim13

echo $SLURM_ARRAY_TASK_ID

Rscript --vanilla SimMain_Script.R $SLURM_ARRAY_TASK_ID

cp -r ./*.RDS "$HOME"/mlVARGD_Sim13


