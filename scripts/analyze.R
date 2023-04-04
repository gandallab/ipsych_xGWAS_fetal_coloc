#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Analyze eCAVIAR results")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gene", help = "")
p <- add_argument(p, "--annot", help = "tri xQTL")


args <- parse_args(p)

# read locus table
test_table <- read.table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv", header = T)
snp <- test_table[args$locus, 'SNP']
chr <- test_table[args$locus, 'CHR']
bp <- test_table[args$locus, 'BP']
gwas <- test_table[args$locus, 'GWAS']

# read ecaviar results
if (grepl("tri", args$annot, fixed = TRUE)) {
    res <- fread(paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_ecaviar_", args$annot, "_col"
), data.table = F)

if (max(res$CLPP > 0.01)) {
    sig <- res %>% filter(CLPP > 0.01)
    sig$locus <- args$locus
    sig$gene <- args$gene
    sig$gwas <- gwas
    write.table(sig, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_ecaviar_", args$annot, "_col_sig.txt"
    ),
    col.names = T, row.names = F, quote = F, sep = "\t"
    )
}
} else {
    res <- fread(paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_ecaviar_col"
), data.table = F)

if (max(res$CLPP > 0.01)) {
    sig <- res %>% filter(CLPP > 0.01)
    sig$locus <- args$locus
    sig$gene <- args$gene
    sig$gwas <- gwas
    write.table(sig, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_ecaviar_col_sig.txt"
    ),
    col.names = T, row.names = F, quote = F, sep = "\t"
    )
}
}
