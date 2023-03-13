#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Subset locus LD to the tested variants")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gene", help = "")

args <- parse_args(p)

# read locus table
test_table <- read.table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv", header = T)
snp <- test_table[args$locus, 'SNP']
chr <- test_table[args$locus, 'CHR']
bp <- test_table[args$locus, 'BP']
gwas <- test_table[args$locus, 'GWAS']

# read full locus LD matrix
if (bp-1e6 > 0) {
    locus_ld <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_", bp-1e6, "_", bp+1e6, ".ld.gz"), 
    data.table = F
)
} else {
    locus_ld <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_0_", bp+1e6, ".ld.gz"), 
    data.table = F
)
}

if (bp-1e6 > 0) {
    locus_bim <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_", bp-1e6, "_", bp+1e6, ".bim"), 
    data.table = F
)
} else {
    locus_bim <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_0_", bp+1e6, ".bim"), 
    data.table = F
)
}

# read list of common variants with QTL
snps <- read.table(paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_snps.txt"
), header = F)
# make sure same order of variants: QTL/GWAS zscores, snp list, LD matrices
# > test <- snps %>% inner_join(locus_bim, by = c("V1"="V2"))
# > is.unsorted(test$V4)
# [1] FALSE

# subset LD matrix
ld <- locus_ld[locus_bim$V2 %in% snps$V1,locus_bim$V2 %in% snps$V1]
write.table(ld, paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_gwas.ld"
), col.names = F, row.names = F, quote = F, sep = "\t")