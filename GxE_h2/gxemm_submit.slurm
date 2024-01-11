#!/bin/bash

#SBATCH --job-name=free
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
#SBATCH --time=1:30:00
#SBATCH --qos=1wk
#SBATCH --array=1-2762
#SBATCH --output=logs/gxemm_%A_%a.out
#SBATCH --error=logs/gxemm_%A_%a.err

module load R
cd /Genomics/argo/users/damelo/projects/HS-Expression-GxE/GxE_h2
Rscript GxEMM-HS_C.R --gene $SLURM_ARRAY_TASK_ID --tissue head --model free --gene_offset 4000
