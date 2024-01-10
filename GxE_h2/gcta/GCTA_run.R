library(optparse)
library(bracer)
library(stringr)

if (sys.nframe() == 0L) {
    option_list <- list(
                        make_option("--gene",
                                    help = ("Gene number to run"),
                                    metavar = "gene"),
                        make_option("--gene_offset",
                                    default = "0",
                                    help = ("Gene offset for running in cluster"),
                                    metavar = "gene_offset"),
                        make_option("--treatment",
                                    help = ("hs, ctrl, hsctrl"),
                                    default = "hs",
                                    metavar = "treatment"),
                        make_option("--outputDir",
                                    default = "/Genomics/ayroleslab2/diogro/projects/NEX-HS_C-GxE/data/output/gcta/",
                                    help = ("OutPut directory"),
                                    metavar = "outputDir")
    )
    parser_object <- OptionParser(usage = "Rscript %prog [Options] --gene [gene]\n",
                                  option_list = option_list,
                                  description = "Run GxEMM")

    opt <- parse_args(parser_object, args = commandArgs(trailingOnly = TRUE),
                      positional_arguments = TRUE)
    i_gene <- as.numeric(opt$options$gene)
    gene_offset <- as.numeric(opt$options$gene_offset)
    treatment <- opt$options$treatment
    out.dir <- opt$options$outputDir

    # quit on error when run non-interactively #small change because this is killing my local sessions T_T
    if (!interactive()) options(error = function() quit(save = "no", status = 1))
}

i_gene = i_gene + gene_offset
print(i_gene)

gene_list = read.csv("gene_list.csv", header = FALSE)


out_name =  paste0(str_pad(i_gene, 4, pad = "0"), "_", gene_list[i_gene, 1], "_", treatment)
out_file = file.path(out.dir, paste0(out_name, ".hsq"))
if(file.exists(out_file)) {
    print("File exists.")
    stop("Out file exists.")
}

if(!file.exists(out.dir)){
    dir.create(out.dir, showWarnings = TRUE, recursive = TRUE)
}
out_file = file.path(out.dir, out_name)
if(treatment == "hs" | treatment == "ctrl"){
    command  = paste0("gcta --reml --grm grm --mpheno ", i_gene,
                      " --pheno VOOMCounts_CPM1_head_", treatment,
                      "_covfree_4svs_CORRECT_Jan8.21.txt --out ", out_file)  
} else if(treatment == "hsctrl"){
    command  = paste0("gcta --reml --grm grm --mpheno ", i_gene,
                      " --gxe HSC.gxe",
                      " --pheno VOOMCounts_CPM1_head_", treatment,
                      "_covfree_4svs_CORRECT_Jan8.21.txt --out ", out_file)  
} else
    stop("Unknown treatment. Use: hs, ctrl or hsctrl for a GxE model.")
print(command)
system(command)
