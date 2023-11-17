packages <- c("gdsfmt", 
              "SeqArray", 
              "SNPRelate", 
              "rio", 
              "plyr",
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
              "viridis",
              "patchwork", 
              "ggrepel", 
              "clusterProfiler", 
              "pointblank",
              "apeglm",
              "ashr", "dcgerard/vicar",
              "furrr",
              "qvalue",
              "ggblend")
# pak::pkg_install(packages)
packages[grepl("/", packages)] <- gsub(".+/", "", packages[grepl("/", packages)])
lapply(packages, library, character.only = TRUE)


