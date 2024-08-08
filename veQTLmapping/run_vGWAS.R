source(here::here("Rscripts/functions.R"))

future::plan(
  list(future::tweak(future::multisession, workers = 32))
  )

if(!require('vGWAS')) {
  install.packages('vGWAS', repos = 'http://r-forge.r-project.org')
  library('vGWAS')
}

rnaseq = import(here::here("cache/rnaseq_all_2024-03-21.rds"), trust = TRUE)

t = 'head'

phenos = t(rnaseq[[t]]$l2c)
rownames(phenos)
covariates = rnaseq[[t]]$covariates

Xp = read.plink(paste0("eQTLmapping/bed_files/", t))
X = as(Xp$genotypes,'numeric') 
rownames(X) = Xp$fam$member
colnames(X) = Xp$map$snp.name

phenos = phenos[rownames(phenos) %in% rownames(X),]

stopifnot(all(rownames(phenos) %in% rownames(X)))
stopifnot(all(rownames(X) %in% rownames(phenos)))

colnames(rnaseq[[t]]$mwash$l2c$Zhat) = paste0("Zhat", 1:dim(rnaseq[[t]]$mwash$l2c$Zhat)[2])
input_df = cbind(covariates, 
                 rnaseq[[t]]$mwash$l2c$Zhat)

input_df <- input_df[match(rownames(phenos), input_df$id),]

stopifnot(rownames(phenos) == input_df$id)

data(pheno)
data(geno)
data(chr)
data(map)

vgwa <- vGWAS(phenotype = pheno, geno.matrix = geno, marker.map = map, chr.index = chr)

vGWAS.variance(phenotype = pheno,
 marker.genotype = geno[,vgwa$p.value == min(vgwa$p.value)])
