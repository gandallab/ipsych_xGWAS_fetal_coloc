#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Get eGene in locus")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--annot", help = "eqtl/isoqtl/sqtl/tri")

args <- parse_args(p)

test_table <- read.table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv", header = T)

if (args$annot == "eqtl") {
    egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/significant_assoc.txt",
    data.table = F
)
} else if (args$annot == "isoqtl") {
    egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/mixed_nominal_70hcp/significant_assoc.txt",
    data.table = F
)
} else if (args$annot == "sqtl") {
    egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/significant_assoc.txt",
    data.table = F
)
} else if (args$annot == "eqtl_tri1") {
    egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T1-nominal_significant_assoc.txt",
    data.table = F
)
} else if (args$annot == "eqtl_tri2") {
    egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T2-nominal_significant_assoc.txt",
    data.table = F
)
} else if (args$annot == "isoqtl_tri1") {
    egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri1_nominal_35hcp/significant_assoc.txt",
    data.table = F
)
} else if (args$annot == "isoqtl_tri2") {
    egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri2_nominal_20hcp/significant_assoc.txt",
    data.table = F
)
} else if (args$annot == "sqtl_tri1") {
    egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri1_nominal_15hcp/significant_assoc.txt",
    data.table = F
)
} else if (args$annot == "sqtl_tri2") {
    egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri2_nominal_10hcp/significant_assoc.txt",
    data.table = F
)
}



snp <- test_table[args$locus, 'SNP']
chr <- test_table[args$locus, 'CHR']
bp <- test_table[args$locus, 'BP']
gwas <- test_table[args$locus, 'GWAS']

if (bp-1e6 > 0) {
    bim <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_", bp-1e6, "_", bp+1e6, ".bim"), 
    data.table = F
)
} else {
    bim <- fread(
    paste0("/u/project/gandalm/cindywen/ipsych_gwas/data/index_snps_ld_matrices/IndexSNPsRegions_", chr, "_0_", bp+1e6, ".bim"), 
    data.table = F
)
}

egene_in_locus <- egene %>% filter(sid %in% bim$V2)

if (length(unique(egene_in_locus$pid)) > 0) {
    egene_list <- data.frame(
    'locus' = args$locus, 
    'gene' = unique(egene_in_locus$pid)
    )

    if (args$annot == "eqtl") {
        write.table(egene_list, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_egene.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "isoqtl") {
        write.table(egene_list, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_isoform.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "sqtl") {
        write.table(egene_list, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_intron.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "eqtl_tri1") {
        write.table(egene_list, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_egene_tri1.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "eqtl_tri2") {
        write.table(egene_list, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_egene_tri2.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "isoqtl_tri1") {
        write.table(egene_list, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_isoform_tri1_35hcp.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "isoqtl_tri2") {
        write.table(egene_list, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_isoform_tri2_20hcp.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "sqtl_tri1") {
        write.table(egene_list, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_intron_tri1_15hcp.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "sqtl_tri2") {
        write.table(egene_list, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_intron_tri2_10hcp.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    }
} else {
    x <- data.frame()
    if (args$annot == "eqtl") {
        write.table(x, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_egene.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "isoqtl") {
        write.table(x, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_isoform.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "sqtl") {
        write.table(x, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_intron.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "eqtl_tri1") {
        write.table(x, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_egene_tri1.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "eqtl_tri2") {
        write.table(x, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_egene_tri2.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "isoqtl_tri1") {
        write.table(x, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_isoform_tri1_35hcp.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "isoqtl_tri2") {
        write.table(x, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_isoform_tri2_20hcp.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "sqtl_tri1") {
        write.table(x, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_intron_tri1_15hcp.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    } else if (args$annot == "sqtl_tri2") {
        write.table(x, paste0(
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus", args$locus, "/locus_intron_tri2_10hcp.txt"
    ), col.names = F, row.names = F, quote = F, sep = "\t")
    }
}

