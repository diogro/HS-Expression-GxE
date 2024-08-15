#!/bin/bash

#sbatch --job-name=gcta
#sbatch --mem=4G --cpus-per-task=1
#sbatch --time=1:30:00
#sbatch --qos=1wk
#SBATCH --array=1-1000
#SBATCH --output=out/gcta_%A_%a.out
#SBATCH --error=out/gcta_%A_%a.err

cd /Genomics/argo/users/damelo/projects/HS-Expression-GxE/gcta
Rscript GCTA_run.R --gene $SLURM_ARRAY_TASK_ID --treatment ctrl --gene_offset 6000
Rscript GCTA_run.R --gene $SLURM_ARRAY_TASK_ID --treatment hs --gene_offset 6000
