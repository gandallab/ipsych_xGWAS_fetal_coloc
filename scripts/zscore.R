#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Extract overlapping variants, write zscore files")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gene", help = "")

args <- parse_args(p)

# read locus table
test_table <- read.table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv", header = T)
snp <- test_table[args$locus, 'SNP']
chr <- test_table[args$locus, 'CHR']
bp <- test_table[args$locus, 'BP']
gwas <- test_table[args$locus, 'GWAS']


# read gene cis assoc, gwas sumstats, locus variant BIM
qtl <- fread(paste0("/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_all_pairs.txt"), data.table = F)
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
gwas_in_qtl <- gwas_sumstats %>% filter(SNP %in% qtl$V2)
qtl_in_gwas <- qtl_in_gwas %>%
    mutate(zscore = V8 / V9)
# > is.unsorted(gwas_in_qtl$BP)
# [1] FALSE
# > sum(qtl_in_gwas$V2 == gwas_in_qtl$SNP)
# [1] 1645

# # check that effect allele is the same between QTL and GWAS
# qtl_bim <- fread("/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test.bim", data.table=F)
# qtl_in_gwas <- qtl_in_gwas %>% inner_join(qtl_bim, by = c("V2"="V2"))
# > is.unsorted(qtl_in_gwas$V4.y)
# [1] FALSE
# > sum(qtl_in_gwas$V2 == gwas_in_qtl$SNP)
# [1] 1645

qtl_in_gwas <- qtl_in_gwas %>% select(V2, zscore)
gwas_in_qtl <- gwas_in_qtl %>% select(SNP, Z)
write.table(qtl_in_gwas, paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_zscore.txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)

write.table(gwas_in_qtl, paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_gwas_zscore.txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)

write.table(as.data.frame(gwas_in_qtl$SNP), paste0(
    "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/", args$gene, "_snps.txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)