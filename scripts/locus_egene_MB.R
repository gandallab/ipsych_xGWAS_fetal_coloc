#! /usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Get eGene in locus")
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

# 1. read GWAS locus bim and LD file
if (bp_19-1e6 > 0) {
    bim <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_", bp_19-1e6, "_", bp_19+1e6, ".bim"), 
    data.table = F
)
    locus_ld <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_", bp_19-1e6, "_", bp_19+1e6, ".ld.gz"), 
    data.table = F
)
} else {
    bim <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_0_", bp_19+1e6, ".bim"), 
    data.table = F
)
    locus_ld <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_0_", bp_19+1e6, ".ld.gz"), 
    data.table = F
)
}

# 2. read GWAS assoc; check alleles between bim (LD) and GWAS, flip Z if neccessary
# A1 GWAS effect; A2 GWAS other
# V5 BIM minor (LD effect allele); V6 BIM major
gwas_sumstats <- fread(paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/sumstats_filtered/iPSYCH2015_EUR_", gwas, ".assoc"), data.table = F)
gwas_sumstats <- gwas_sumstats %>% inner_join(bim, by = c("SNP" = "V2"))
gwas_sumstats <- gwas_sumstats %>% filter((A1 == V5 & A2 == V6) | (A1 == V6 & A2 == V5))
gwas_sumstats[which(gwas_sumstats$A1 == gwas_sumstats$V6 & gwas_sumstats$A2 == gwas_sumstats$V5),'Z'] <- -gwas_sumstats[which(gwas_sumstats$A1 == gwas_sumstats$V6 & gwas_sumstats$A2 == gwas_sumstats$V5),'Z']


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

# 6. read 1KG EUR variants, need to restrict to these to calculate LD for MB QTL
eur_1kg <- fread(paste0("/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.", chr, ".bim"), data.table = F)



# 7. write ecaviar input files for candidate genes: zscore for QTL and GWAS for shared variants; LD from GWAS
if (length(unique(sig_locus$Gene)) > 0) {
    egene_list <- data.frame(
    'locus' = args$locus, 
    'gene' = unique(sig_locus$Gene)
    )
    
    for(gene in unique(sig_locus$Gene)) {
        gene_full_assoc <- full_assoc %>% filter(Gene == gene)
        gene_full_assoc <- gene_full_assoc %>% 
        separate(SNP.x, c("snpchr", "snppos", "rsid", "snpalleles"), sep = ":") %>%
        separate(SNPAlleles.x, c("A1", "A2"), sep = "/")
        # MB has multiple rows with same rsid and multiple alleles, check against GWAS and keep the entry with matching alleles
        shared <- gene_full_assoc %>% inner_join(gwas_sumstats, by = c("rsid" = "SNP"))
        # A1.x: QTL other allele; A2.x: QTL effect allele
        # A1.y: GWAS effect allele; A2.y: GWAS other allele
        # V5: GWAS bim minor; V6: GWAS bim major
        # GWAS z respective to V5/V6 if inconsistent between GWAS sum stats and bim
        shared <- shared %>% filter((A1.x == A2.y & A2.x == A1.y) | (A1.x == A1.y & A2.x == A2.y))

        # order by GWAS hg19
        shared <- shared %>% arrange(BP) 

        # restrict to 1KG EUR variants, and check alleles
        # V5.y: EUR 1KG minor; V6.y: EUR 1KG major
        shared <- shared %>% inner_join(eur_1kg, by = c("rsid" = "V2"))
        shared <- shared %>% filter((A1.x == V6.y & A2.x == V5.y) | (A1.x == V5.y & A2.x == V6.y))
        shared[which(shared$A1.x == shared$V5.y & shared$A2.x == shared$V6.y),'MetaPZ.x'] <- -shared[which(shared$A1.x == shared$V5.y & shared$A2.x == shared$V6.y),'MetaPZ.x']


        # eCAVIAR input files
        qtl_z <- shared %>% select(rsid, MetaPZ.x)
        gwas_z <- shared %>% select(rsid, Z)
        snps <- as.data.frame(shared$rsid)
        colnames(snps)[1] <- "V1"
        ld <- locus_ld[bim$V2 %in% snps$V1, bim$V2 %in% snps$V1]
        
        write.table(qtl_z, paste0(gene, "_MB_zscore.txt"), col.names = F, row.names = F, quote = F, sep = "\t")

        write.table(gwas_z, paste0(gene, "_MB_gwas_zscore.txt"), col.names = F, row.names = F, quote = F, sep = "\t")

        write.table(snps, paste0(gene, "_MB_snps.txt"), col.names = F, row.names = F, quote = F, sep = "\t")

        write.table(ld, paste0(gene, "_MB_gwas.ld"), col.names = F, row.names = F, quote = F, sep = "\t")
    }
} else {
    egene_list <- data.frame()
}


write.table(egene_list, paste0("locus_egene_MB.txt"), col.names = F, row.names = F, quote = F, sep = "\t")