#!/bin/bash

#SBATCH --job-name=gxemm
#SBATCH --mem=8G 
#SBATCH --cpus-per-task=1
#SBATCH --time=1:30:00
#SBATCH --qos=1wk
#SBATCH --array=3-5
#SBATCH --output=logs/gxemm_%A_%a.out
#SBATCH --error=logs/gxemm_%A_%a.err

module load R
cd /Genomics/argo/users/damelo/projects/HS-Expression-GxE/GxE_h2
Rscript GxEMM-HS_C.R --gene $SLURM_ARRAY_TASK_ID --tissue head --model iid 
