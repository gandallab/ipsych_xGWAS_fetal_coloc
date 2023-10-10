from os.path import join
import os
import numpy as np
import pandas as pd
import sys

TEST_TABLE = pd.read_table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP_hg38.tsv")
LOCUS_GENE_TABLE = pd.read_table("locus_gene_list_MB.txt")

def get_1kg_eur_bfile_chr(wildcards):
    TEST_TABLE = pd.read_table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP_hg38.tsv").set_index(
        "locus", drop=True
    )
    chromosome = TEST_TABLE.loc[int(wildcards.locus), "CHR"]
    return (
        "/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC."
        + str(chromosome)
    )

"""
rules:
- locus_egene: for each GWAS locus, get a list of MB gene to run eCAVIAR; also write eCAVIAR input files for candidate genes: zscores, LD
- write_locus_gene_list: concat list of candidate genes from all GWAS loci
- ld_qtl: calculate LD matrix for MB eQTL
- ecaviar
- analyze_ecaviar
- concat
- coloc
"""
rule all:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/MB_CLPP_sig.txt",
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/coloc_MB.done",
            locus=TEST_TABLE.locus.values,
        ),
        
rule locus_egene:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP_hg38.tsv",
        "/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-TopEffects.txt.gz",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_egene_MB.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/locus_egene_MB.R \
            --locus {wildcards.locus}
        """

rule write_locus_gene_list:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_egene_MB.txt",
            locus=TEST_TABLE.locus.values,
        ),
    output:
        "locus_gene_list_MB.txt",
    shell:
        """
        cat /u/project/gandalm/cindywen/ipsych_gwas/out/locus*/locus_egene_MB.txt > {output[0]}
        sed  -i '1i locus\tgene' {output[0]}
        """

rule ld_qtl:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB.ld",
    params:
        get_1kg_eur_bfile_chr,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile {params[0]} \
            --r \
            --matrix \
            --extract {input[0]} \
            --out /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}_MB
        """

rule ecaviar:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB_gwas.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB_gwas_zscore.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB_ecaviar_col",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_MB_ecaviar \
                -l {wildcards.gene}_MB_gwas.ld \
                -z {wildcards.gene}_MB_gwas_zscore.txt \
                -l {wildcards.gene}_MB.ld \
                -z {wildcards.gene}_MB_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB_ecaviar_col",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB.done"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/analyze.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot MB
        touch {output[0]}
        """

rule concat:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_MB.done",
            zip,
            locus = LOCUS_GENE_TABLE.locus.values,
            gene = LOCUS_GENE_TABLE.gene.values,
        ),
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/MB_CLPP_sig.txt",
    shell:
        """
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' locus*/*MB_ecaviar_col_sig.txt > {output[0]}
        """

rule coloc:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP_hg38.tsv",
        "/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-TopEffects.txt.gz",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/coloc_MB.done",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/coloc_MB.R \
            --locus {wildcards.locus}
        touch {output[0]}
        """