#! /usr/bin/env Rscript

library(coloc)
library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Run coloc")
p <- add_argument(p, "--locus", help = "")
args <- parse_args(p)


test_table <- read.table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP_hg38.tsv", header = T)

snp <- test_table[args$locus, 'SNP']
# chr_38 <- test_table[args$locus, 'chr_38']
bp_38 <- test_table[args$locus, 'bp_38']
chr <- test_table[args$locus, 'CHR']
bp_19 <- test_table[args$locus, 'BP']
gwas <- test_table[args$locus, 'GWAS']


setwd(paste0("/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/"))

# 1. read gwas variants
if (bp_19-1e6 > 0) {
    bim <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_", bp_19-1e6, "_", bp_19+1e6, ".bim"), 
    data.table = F
)
#     locus_ld <- fread(
#     paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_", bp_19-1e6, "_", bp_19+1e6, ".ld.gz"), 
#     data.table = F
# )
} else {
    bim <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_0_", bp_19+1e6, ".bim"), 
    data.table = F
)
#     locus_ld <- fread(
#     paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_0_", bp_19+1e6, ".ld.gz"), 
#     data.table = F
# )
}

# 2. read GWAS assoc
gwas_sumstats <- fread(paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/sumstats_filtered/iPSYCH2015_EUR_", gwas, ".assoc"), data.table = F)
gwas_sumstats <- gwas_sumstats %>% inner_join(bim, by = c("SNP" = "V2"))

# 3. read MB full assoc file
full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-chr", chr, ".txt.gz"), data.table = F, sep = "\t")
perm <- fread("/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-TopEffects.txt.gz", data.table = F, sep = "\t")
full_assoc <- full_assoc %>% select(-c(18:20))
perm <- perm %>% select(-c(18:20))
full_assoc <- full_assoc %>% left_join(perm, by = "Gene")

# 4. filter for nominal significant eQTL
sig_assoc <- full_assoc %>% filter(MetaP.x < PvalueNominalThreshold)
sig_assoc <- sig_assoc %>% separate(SNP.x, c("snpchr", "snppos", "rsid", "snpalleles"), sep = ":")

# 5. find genes with GWAS variants as eQTL
sig_locus <- sig_assoc %>% filter(rsid %in% gwas_sumstats$SNP)

# 6. run coloc. Flip beta where neccessary
if (length(unique(sig_locus$Gene)) > 0) {
    for(gene in unique(sig_locus$Gene)) {
        gene_full_assoc <- full_assoc %>% filter(Gene == gene)
        gene_full_assoc <- gene_full_assoc %>% 
        separate(SNP.x, c("snpchr", "snppos", "rsid", "snpalleles"), sep = ":") %>%
        separate(SNPAlleles.x, c("A1", "A2"), sep = "/")
        shared <- gene_full_assoc %>% inner_join(gwas_sumstats, by = c("rsid" = "SNP"))
        shared <- shared %>% filter((A1.x == A2.y & A2.x == A1.y) | (A1.x == A1.y & A2.x == A2.y))
        # flipped alleles, flip QTL zscore
        shared[which(shared$A1.x == shared$A1.y & shared$A2.x == shared$A2.y),'MetaBeta.x'] <- -shared[which(shared$A1.x == shared$A1.y & shared$A2.x == shared$A2.y),'MetaBeta.x']
        # creat coloc data
    feature_data <- list("beta" = shared$MetaBeta.x,
                     "varbeta" = shared$MetaSE.x * shared$MetaSE.x,
                     "snp" = shared$rsid,
                     "position" = shared$BP,
                     "type" = "quant",
                     "sdY" =1)                                   
    gwas_data <- list("beta" = shared$BETA,
                  "varbeta" = shared$SE * shared$SE,
                  "snp" = shared$rsid,
                  "position" = shared$BP,
                  "type" = "quant",
                  "sdY" = 1)
    res <- coloc.abf(dataset1 = feature_data, dataset2 = gwas_data)
    saveRDS(res, paste0(gene, ".MB.coloc.res.rds"))
    if (res$summary[6] > 0.7) {
        write.table(res$summary, paste0(gene, ".MB.coloc.sigPP4"))
    }

    }
} else {
    # exit gracefully if no candidate target features
    quit(save = "no", status = 0)
}
