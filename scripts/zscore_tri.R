#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Extract overlapping variants, write zscore files. Added annot for file extensions; also need to calculated zscore from beta and pval for trimester QTL")
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


# read feature cis assoc
# add fill = T as some rows have missing FDR column in PEC eQTL file, fread terminate
qtl <- fread(paste0("/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_all_pairs_", args$annot, ".txt"), data.table = F, fill = T)

# read locus GWAS
gwas_sumstats <- fread(paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/sumstats_filtered/iPSYCH2015_EUR_", gwas, ".assoc"), data.table = F)
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

# subset GWAS to locus
gwas_sumstats <- gwas_sumstats %>% filter(SNP %in% locus_bim$V2)
# > is.unsorted(gwas_sumstats$BP)
# [1] TRUE

# order GWAS variants as in BIM
gwas_sumstats <- gwas_sumstats[order(match(gwas_sumstats[, 'SNP'], locus_bim[, 'V2'])), ]
# > is.unsorted(gwas_sumstats$BP)
# [1] FALSE

# overlap variants QTL and GWAS locus
qtl_in_gwas <- qtl %>% filter(V2 %in% gwas_sumstats$SNP)
# somehow tri sQTL has duplicate variants; remove these, all clu_38813 introns
qtl_in_gwas <- qtl_in_gwas[!duplicated(qtl_in_gwas$V2),]
gwas_in_qtl <- gwas_sumstats %>% filter(SNP %in% qtl$V2)
# order QTL. Note isoQTL and sQTL from GTEx's FastQTL variants are not ordered as in BIM. Also comfirmed same BP between data_bim and locus_bim
# qtl_in_gwas <- qtl_in_gwas %>% inner_join(locus_bim, by = c("V2" ="V2")) # this is unneccessary
qtl_in_gwas <- qtl_in_gwas[order(match(qtl_in_gwas[, 'V2'], locus_bim[, 'V2'])), ]


qtl_in_gwas <- qtl_in_gwas %>%
    mutate(zscore = sign(V5) * abs(qnorm(V4 / 2))) %>%
    select(V2, zscore)



qtl_in_gwas <- qtl_in_gwas %>% select(V2, zscore)
gwas_in_qtl <- gwas_in_qtl %>% select(SNP, Z)


write.table(qtl_in_gwas, paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_zscore_", args$annot, ".txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)

write.table(gwas_in_qtl, paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_gwas_zscore_", args$annot, ".txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)

write.table(as.data.frame(gwas_in_qtl$SNP), paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_snps_", args$annot, ".txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)

