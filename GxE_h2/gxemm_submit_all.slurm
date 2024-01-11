#!/bin/bash

#SBATCH --job-name=gxe-all
#SBATCH --mem=12G 
#SBATCH --cpus-per-task=3
#SBATCH --time=1:30:00
#SBATCH --qos=1wk
#SBATCH --array=1-2448
#SBATCH --output=logs/gxemm_%A_%a.out
#SBATCH --error=logs/gxemm_%A_%a.err

module load R
cd /Genomics/argo/users/damelo/projects/HS-Expression-GxE/GxE_h2
Rscript GxEMM-HS_C.R --gene $SLURM_ARRAY_TASK_ID --tissue body --model hom  --gene_offset 4000 &
Rscript GxEMM-HS_C.R --gene $SLURM_ARRAY_TASK_ID --tissue body --model iid  --gene_offset 4000 &
Rscript GxEMM-HS_C.R --gene $SLURM_ARRAY_TASK_ID --tissue body --model free --gene_offset 4000
