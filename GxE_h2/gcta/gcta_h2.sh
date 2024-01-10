#!/bin/bash

gcta --reml \
    --grm grm \
    --mpheno 2 \
    --pheno /Genomics/grid/users/damelo/projects/NEX-HS_C-GxE/gcta/VOOMCounts_CPM1_head_ctrl_covfree_4svs_CORRECT_Jan8.21.txt \
    --out test_pheno
