#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Analyze eCAVIAR results")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gene", help = "")
p <- add_argument(p, "--annot", help = "")


args <- parse_args(p)

if (args$annot == "MB") {
    # MetaBrain eQTL
    test_table <- read.table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP_hg38.tsv", header = T)
    snp <- test_table[args$locus, 'SNP']
    # chr_38 <- test_table[args$locus, 'chr_38']
    bp_38 <- test_table[args$locus, 'bp_38']
    chr <- test_table[args$locus, 'CHR']
    bp_19 <- test_table[args$locus, 'BP']
    gwas <- test_table[args$locus, 'GWAS']
    res <- fread(paste0("/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_MB_ecaviar_col"), data.table = F)
    if (max(res$CLPP > 0.01)) {
        sig <- res %>% filter(CLPP > 0.01)
        sig$locus <- args$locus
        sig$gene <- args$gene
        sig$gwas <- gwas
        write.table(sig, paste0("/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_MB_ecaviar_col_sig.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
}
} else {

test_table <- read.table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv", header = T)
snp <- test_table[args$locus, 'SNP']
chr <- test_table[args$locus, 'CHR']
bp <- test_table[args$locus, 'BP']
gwas <- test_table[args$locus, 'GWAS']

if (grepl("(tri|pec)", args$annot)) {
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
} else if (args$annot == "thistle") {
    # thistle sQTL
    res <- fread(paste0("/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_thistle_ecaviar_col"), data.table = F)
    if (max(res$CLPP > 0.01)) {
    sig <- res %>% filter(CLPP > 0.01)
    sig$locus <- args$locus
    sig$gene <- args$gene
    sig$gwas <- gwas
    write.table(sig, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_thistle_ecaviar_col_sig.txt"
    ),
    col.names = T, row.names = F, quote = F, sep = "\t"
    )
}

} else {
    # fetal e/iso/sQTL
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
}