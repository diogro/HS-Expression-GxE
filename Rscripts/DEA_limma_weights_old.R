
#conda activate limma
#check that R version is 4.


library(limma)
library(edgeR) 

start=Sys.time()
Sys.time()

cat("loading data \n")
countdata <- read.table("data/GXEpaper/genecounts/RawCounts_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt", h=T, check.names=F)
coldata <- read.table("data/GXEpaper/genecounts/Covariates_forGEMMA_Jan82021.txt",h=T)
coldata$treatment <- as.factor(coldata$treatment)
coldata$RNAlibBatch <- as.factor(coldata$RNAlibBatch)
coldata$RNAseqBatch <- as.factor(coldata$RNAseqBatch)
coldata$egglayBatch <- as.factor(coldata$egglayBatch)
coldata$platingBatch <- as.factor(coldata$platingBatch)
cat ("data loaded \n")
cat("are samples ordered in the right way in coldata and countdata matrices?", all.equal(rownames(coldata),colnames(countdata)), "\n")
cat (nrow(coldata), "samples, ", nrow(countdata), "genes \n")

design <- model.matrix(~0+coldata$treatment+
                         coldata$sv1+coldata$sv2+coldata$sv3+coldata$sv4+
                         coldata$egglayBatch+coldata$RNAseqBatch+coldata$platingBatch+coldata$RNAlibBatch) #ctrl=1, hs=2 
colnames(design) <- c("ctrl","hs", "sv1", "sv2", "sv3","sv4", "e2", "e3", "s2", "s3", "s4", "s5","p22","p23","p24","p25", "p26", "p27", "l2", "lb3")
contrast <- makeContrasts(HSvsC=hs-ctrl, levels = design)

tmp1<-Sys.time()
cat("normalizing and estimating mean-variance weights \n")
countdata.list <- DGEList(counts=countdata)
countdata.norm <- calcNormFactors(countdata.list) #gets norm. factors based on TMM (controlling not only for lib size, but also composition)
countdata.voom <- voom(countdata.norm, design = design, plot=T) #the E matrix of normalized counts is the same as in boomWithQualityWeights
cat("running analysis \n")
fit <- lmFit(countdata.voom,design)
fitc <- contrasts.fit(fit, contrasts = contrast)
fit2 <- eBayes(fitc)
summary(decideTests(fit2))
res<- topTable(fit2, number = Inf)
tmp3<-Sys.time()
cat("analysis done", tmp3-tmp1, "\n")

cat("writing results table and RDS object with all info from fit2 object \n")
write.table(res, "DiffExpression_results_head_HSvsCTRL_March21.txt", col.names = T, row.names = T, quote = F, sep="\t")
saveRDS(fit2, "DiffExpression_results_head_HSvsCTRL_March21.rds")
saveRDS(fit, "DiffExpression_lmFit_head_HSvsCTRL_March21.rds") #from here I get the mean expression logcpm per group fit$coefficients
#and the stderror fit$stdev.unscaled * fit$sigma

Sys.time()
end=Sys.time()
cat("total time = ", end-start)

x <- removeBatchEffect(countdata.voom,
                       model.matrix(~1+coldata$egglayBatch+coldata$RNAseqBatch+coldata$platingBatch+coldata$RNAlibBatch), 
                       design=model.matrix(~0+coldata$treatment))  



y <- DGEList(countdata)
y <- calcNormFactors(y)
design <- model.matrix(~0+coldata$treatment) # 2 groups
v <- voom(y, design)
covariates <- model.matrix(~coldata$egglayBatch+coldata$RNAseqBatch+coldata$platingBatch+coldata$RNAlibBatch) #ctrl=1, 
no.sv1 <- removeBatchEffect(v, design = design, covariates = covariates)

