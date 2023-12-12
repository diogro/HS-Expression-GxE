library(GridLMM)
library(vcfR)

vcf = read.vcfR("eQTLmapping/vcf_files/imputed/head.imputed.vcf.gz")

gt <- extract.gt(vcf)
gt[1:10, 1:10]

X = gt2num(gt)

map = data.frame(snp = rownames(X$genomat), 
                 chr = vcf@fix[,"CHROM"], 
                 pos = vcf@fix[,"POS"])

data = rnaseq_data$head



input_df = data$covariates |>
     mutate(y = t(data$l2c)[,2]) |>
     filter(id %in% colnames(X$genomat))
tX = t(X$genomat)

global_formula = y ~ 1 + egglayBatch + RNAseqBatch +  platingBatch + RNAlibBatch + treatment + (1|id)

gxe_gwas = GridLMM_GWAS(
                        formula = global_formula, # the same error model is used for each marker. It is specified similarly to lmer
                        test_formula = ~1 + treatment,              # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
                        reduced_formula = ~1,                 # This is the null model for each test. ~0 means just the error model. ~1 means the null includes the intercept of the marker, but not anything additional
                        data = as.data.frame(input_df),                          # The dataframe to look for terms from the 3 models
                        weights = NULL,                       # optional observation-specific weights
                        X = tX,                                # The matrix of markers. Note: This can be of dimension n_g x p, where n_g is the number of unique genotypes.
                        X_ID = 'id',                        # The column of the data that identifies each genotype. Each level of data$Geno should match a rowname of X
                        h2_start = NULL,                      # An optional vector of h2s to use as starting values for each random effect in the model. If NULL, will be calculated from the error model using GridLMM_ML
                        h2_step = 0.01,                       # step size per random effect for testing alternate values of h2
                        max_steps = 100,                      # maximum number of steps of size h2_step to take from h2_start
                        X_map = map,                          # Optional. The marker positions.
                        relmat = NULL,                        # A list of Kernel matrices for the random effects. If X_ID (here Geno) is not included in this list, then it is calculated as tcrossprod(Xc)/ncol(Xc) where Xc is the centered (and optionally scaled) X. If any random effects are described in `error_model` but not provided here, the Kernel is assumed to be the identity matrix
                        centerX = TRUE,                       # Should the markers be centered when calculating the GRM (only will be done if needed for calculating the GRM),
                        scaleX = FALSE,                       # Should the markers be scaled to have constant variance when calculating the GRM?
                        fillNAX = FALSE,                      # Should missing marker data be filled in with the mean allele frequency?
                        method = 'REML',                      # Should the best model be selected by REML (if False, will be selected by ML)
                        mc.cores = 1,#my_detectCores(),          # How many cores should be used for parallel processing. Unless X is large, tends to actually be faster with mc.cores = 1
                        verbose = T                           # Should progress be printed to the screen?
)

gxe_gwas$results$p_value_REML.1
png("tmp/gxe_gwas.png")
hist(gxe_gwas$results$p_value_REML.2, breaks = 100)
dev.off()
