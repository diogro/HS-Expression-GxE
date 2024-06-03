# pak::pkg_install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
gene_locations_GR <- transcripts(txdb, columns = "gene_id")
gene_locations_GR$gene_id = unlist(gene_locations_GR$gene_id)
chrs <- c("2L", "2R", "3L", "3R", "4", "X")

# snp_pos = data.frame(rs = snp_locations) |>
#   mutate(pos = rs) |>
#   separate(pos, into = c("chr", "bp"), sep = "_") |>
#   mutate(chr = paste0("chr", as.character(chr))) |>
#   mutate(bp = as.numeric(bp))
# export(snp_pos, file = here::here("cache/snp_pos.rds"))
snp_pos = import(here::here("cache/snp_pos.rds"))

getCisSnps = function(transcript, window = 10000, snps = snp_pos){
  snp_GR <- GRanges(
      seqnames = snps$chr,
      ranges = IRanges(snps$bp, end = snps$bp, 
                      names = snps$rs))
  gene_location = gene_locations_GR[which(gene_locations_GR$gene_id == transcript)] + window
  return(snps$rs[unique(findOverlaps(snp_GR, gene_location)@from)])
}