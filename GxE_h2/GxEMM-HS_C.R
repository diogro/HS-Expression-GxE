# pak::pkg_install(c("optparse", "bracer"))
library(optparse)
library(bracer)

if (sys.nframe() == 0L) {
    option_list <- list(
                        make_option("--gene",
                                    help = ("Gene number to run"),
                                    metavar = "gene"),
                        make_option("--gene_offset",
                                    default = "0",
                                    help = ("Gene offset for running in cluster"),
                                    metavar = "gene_offset"),
                        make_option("--tissue",
                                    default = "head",
                                    help = ("Tissue to run (head or body)"),
                                    metavar = "tissue"),
                        make_option("--model",
                                    help = ("Model to run in GxEMM (free, iid or hom)"),
                                    default = "hom",
                                    metavar = "model"),
                        make_option("--input",
                                    help = ("Path to gene expression rds file."),
                                    default = "/Genomics/argo/users/damelo/projects/HS-Expression-GxE/cache/rnaseq_all_2024-03-21.rds",
                                    metavar = "input"),
                        make_option("--grm", 
                                    default = "/Genomics/argo/users/damelo/projects/HS-Expression-GxE/eQTLmapping/GRMs",
                                    help = ("Path to GRM file"),
                                    metavar = "grm"),
                        make_option("--outputDir",
                                    default = "/Genomics/argo/users/damelo/projects/HS-Expression-GxE/output/GxEMM/",
                                    help = ("OutPut directory"),
                                    metavar = "outputDir")
    )
    parser_object <- OptionParser(usage = "Rscript %prog [Options] --gene [gene]\n",
                                  option_list = option_list,
                                  description = "Run GxEMM")

    ## TO TEST INTERACTIVELY the command-line arguments
    # input <- "--gene 1"
    # command.args <- strsplit(input, " ")[[1]]
    # opt <- parse_args(parser_object, args = command.args, positional_arguments = TRUE)
    ## SKIP opt line below
    ## aliases
    opt <- parse_args(parser_object, args = commandArgs(trailingOnly = TRUE),
                      positional_arguments = TRUE)
    i_gene <- as.numeric(opt$options$gene)
    gene_offset <- as.numeric(opt$options$gene_offset)
    tissue <- opt$options$tissue
    model <- opt$options$model
    input_path <- opt$options$input
    GRM_path <- paste0(opt$options$grm, "/", tissue, ".cXX.txt")
    out.dir <- paste0(opt$options$outputDir, "/", tissue, "/")

    # quit on error when run non-interactively #small change because this is killing my local sessions T_T
    if (!interactive()) options(error = function() quit(save = "no", status = 1))
}

library(GxEMM)
library(stringr)
ldak_loc  <- "/Genomics/argo/users/damelo/.local/bin/ldak5.linux"

i_gene = i_gene + gene_offset
print(i_gene)

sample_list = read.table(paste0("/Genomics/argo/users/damelo/projects/HS-Expression-GxE/eQTLmapping/sample_list/", tissue,".txt"), header = FALSE)[,1]

rnaseq_data = readRDS(input_path)
expr_df = t(rnaseq_data[[tissue]]$mwash_residuals$l2c)[sample_list, ]

out_file = file.path(out.dir, paste0(str_pad(i_gene, 4, pad = "0"), "_", colnames(expr_df)[i_gene], "_", model, ".rds"))
if(file.exists(out_file)) {
    print("File exists.")
    stop("Out file exists.")
}

covariates = rnaseq_data[[tissue]]$covariates
covariates = covariates[match(sample_list, covariates$id),]
stopifnot(all(sample_list == covariates$id))
ID_C = covariates[covariates$treatment == 0, "id"]
ID_HS = covariates[covariates$treatment == 1, "id"]

Z_HS = as.numeric(rownames(expr_df) %in% ID_HS)
Z_C = as.numeric(rownames(expr_df) %in% ID_C)

Z_all = cbind(Z_HS, Z_C)
X     <- Z_all[,-1] 

GRM = as.matrix(read.table(GRM_path))

if(!file.exists(out.dir)){
    dir.create(out.dir, showWarnings = TRUE, recursive = TRUE)
}

if(model == "free"){
    out_HSC = GxEMM(expr_df[,i_gene], X = X, K = GRM, Z = Z_all, 
                        gtype='free', etype='free', ldak_loc=ldak_loc)
} else{ 
    out_HSC = GxEMM(expr_df[,i_gene], X = X, K = GRM, Z = Z_all, 
                         gtype=model, ldak_loc=ldak_loc)
} 

if(identical(out_HSC, "LDAK did not converge")){
    print("LDAK did not converge.")
    stop("LDAK did not converge.")
    model = "FAIL"
}
out_file = file.path(out.dir, paste0(str_pad(i_gene, 4, pad = "0"), "_", colnames(expr_df)[i_gene], "_", model, ".rds"))
saveRDS(out_HSC, file = out_file)
