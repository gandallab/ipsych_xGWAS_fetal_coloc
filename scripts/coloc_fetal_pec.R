#! /usr/bin/env Rscript

library(coloc)
library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Run coloc")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--annot", help = "eqtl, sqtl, isoqtl, tri1_eqtl, ... , pec_eqtl, pec_isoqtl")

args <- parse_args(p)


test_table <- read.table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv", header = T)

snp <- test_table[args$locus, 'SNP']
# chr_38 <- test_table[args$locus, 'chr_38']
# bp_38 <- test_table[args$locus, 'bp_38']
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

if (args$annot == "eqtl") {
    # 3. read QTL full assoc file
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    # 4. filter for nominal significant QTL
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "isoqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/mixed_nominal_70hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/mixed_nominal_70hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "sqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/gtex_chunk", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/significant_assoc.txt", data.table = F)
} else if (args$annot == "tri1_eqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/tri1_25_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T1-nominal_significant_assoc.txt", data.table = F)
} else if (args$annot == "tri2_eqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/tri2_15_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T2-nominal_significant_assoc.txt", data.table = F)
} else if (args$annot == "tri1_isoqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri1_nominal_35hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri1_nominal_35hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "tri2_isoqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri2_nominal_20hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri2_nominal_20hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "tri1_sqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri1_nominal_15hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri1_nominal_15hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "tri2_sqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri2_nominal_10hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri2_nominal_10hcp/significant_assoc.txt", data.table = F)
}

colnames(full_assoc) <- c("gene_id", "variant_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se")
# trimester files has NA
full_assoc <- full_assoc[complete.cases(full_assoc),]

# get A1/A2 of QTL sum stats
if (args$annot %in% c("eqtl", "isoqtl", "sqtl")) {
    qtl_bim <- fread("/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test.bim", data.table = F)
} else if (grepl("tri", args$annot)) {
    qtl_bim <- fread("/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeRel.bim", data.table = F)
}
full_assoc <- full_assoc %>% inner_join(qtl_bim, by = c("variant_id" = "V2"))

full_assoc <- full_assoc %>% unite("pair", gene_id, variant_id, sep = "-", remove = FALSE)
full_sig <- full_sig %>% unite("pair", pid, sid, sep = "-", remove = FALSE)
sig_assoc <- full_assoc %>% filter(pair %in% full_sig$pair)
    
# 5. find genes with GWAS variants as QTL
sig_locus <- sig_assoc %>% filter(variant_id %in% gwas_sumstats$SNP)

# 6. run coloc. Flip beta where neccessary
if (length(unique(sig_locus$gene_id)) > 0) {
    for(gene in unique(sig_locus$gene_id)) {
        gene_full_assoc <- full_assoc %>% filter(gene_id == gene)
        shared <- gene_full_assoc %>% inner_join(gwas_sumstats, by = c("variant_id" = "SNP"))
        # QTL and GWAS A1/A2
        shared <- shared %>% filter((V5.x == A1 & V6.x == A2) | (V5.x == A2 & V6.x == A1))
        # flipped alleles, flip QTL beta
        shared[which(shared$V5.x == shared$A2 & shared$V6.x == shared$A1),'slope'] <- -shared[which(shared$V5.x == shared$A2 & shared$V6.x == shared$A1),'slope']
        # creat coloc data
        feature_data <- list("beta" = shared$slope,
                     "varbeta" = shared$slope_se * shared$slope_se,
                     "snp" = shared$variant_id,
                     "position" = shared$BP,
                     "type" = "quant",
                     "sdY" = 1)
        gwas_data <- list("beta" = shared$BETA,
                  "varbeta" = shared$SE * shared$SE,
                  "snp" = shared$variant_id,
                  "position" = shared$BP,
                  "type" = "quant",
                  "sdY" = 1)
    res <- coloc.abf(dataset1 = feature_data, dataset2 = gwas_data)
    saveRDS(res, paste0(gene, ".", args$annot, ".coloc.res.rds"))
    if (res$summary[6] > 0.7) {
        write.table(res$summary, paste0(gene, ".", args$annot, ".coloc.sigPP4"))
    }

    }
} else {
    # exit gracefully if no candidate target features
    quit(save = "no", status = 0)
}