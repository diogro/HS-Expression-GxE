#!/bin/bash

#sbatch --job-name=gcta
#sbatch --mem=4G --cpus-per-task=1
#sbatch --time=1:30:00
#sbatch --qos=1wk
#SBATCH --array=1-8786
#SBATCH --output=out/gcta_%A_%a.out
#SBATCH --error=out/gcta_%A_%a.err

module load R/4.0.3
cd /Genomics/ayroleslab2/diogro/projects/NEX-HS_C-GxE/gcta
Rscript GCTA_run.R --gene $SLURM_ARRAY_TASK_ID --treatment hsctrl --gene_offset 7999
