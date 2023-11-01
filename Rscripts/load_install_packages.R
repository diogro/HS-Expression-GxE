packages <- c("gdsfmt", 
              "SeqArray", 
              "SNPRelate", 
              "rio", 
              "tidyverse", 
              "org.Dm.eg.db", 
              "txtplot",
              "edgeR",
              "limma", 
              "sva", 
              "tictoc",
              "DESeq2",
              "rrcov", 
              "cowplot", 
              "viridis")
# pak::pkg_install(packages)
lapply(packages, library, character.only = TRUE)
