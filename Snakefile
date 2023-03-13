from os.path import join
import os
import numpy as np
import pandas as pd
import sys


# configfile: "config.yaml"

"""
rules:
    - locus_egene: get a list of candidate genes for colocalization
    - write_locus_gene_list: write a list of all locus-gene pairs for parallelization
    - extract_cis_assoc_gene: for each gene, extract all cis assoc in fetal data
    - zscore: for each gene, write zscore files for QTL and GWAS for the common variants; write list of common variants for LD calculation
    - ld_qtl: write LD matrix using fetal QTL genotype data
    - ld_gwas: subset locus LD matrix to common variants
    - ecaviar: run ecaviar
    - analyze: extract CLPP > 0.01
"""

TEST_TABLE = pd.read_table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv")

LOCUS_GENE_TABLE = pd.read_table("locus_gene_list.txt")

# not ideal workaround to avoid wildcards ambiguity
ruleorder: ld_gwas > ld_qtl

rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.done",
            zip,
            locus = LOCUS_GENE_TABLE.locus.values,
            gene = LOCUS_GENE_TABLE.gene.values,
        ),


rule locus_egene:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_egene.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/
        Rscript scripts/locus_egene.R \
            --locus {wildcards.locus}
        """

rule write_locus_gene_list:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_egene.txt",
            locus=TEST_TABLE.locus.values,
        )
    output:
        "locus_gene_list.txt",
    shell:
        """
        cat /u/project/gandalm/cindywen/ipsych_gwas/out/locus*/locus_egene.txt > locus_gene_list.txt
        sed  -i '1i locus\tgene' locus_gene_list.txt
        """

rule extract_cis_assoc_gene:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/gtex.allpairs.txt.gz",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs.txt",
    shell:
        """
        awk -v a="{wildcards.gene}" '$1 == a {{print}}' <(zcat {input[0]}) > {output[0]}
        """

rule zscore:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/zscore.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene}
        """

rule ld_qtl:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test \
            --r \
            --matrix \
            --extract /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}_snps.txt \
            --out /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}
        """

rule ld_gwas:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/ld_gwas.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene}
        """

rule ecaviar:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_col"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_ecaviar \
                -l {wildcards.gene}_gwas.ld \
                -z {wildcards.gene}_gwas_zscore.txt \
                -l {wildcards.gene}.ld \
                -z {wildcards.gene}_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_col",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.done"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/analyze.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene}
        touch {output[0]}
        """